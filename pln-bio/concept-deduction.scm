(define-module (pln-bio concept-deduction)
    #:use-module (opencog)
    #:use-module (opencog logger)
    #:use-module (opencog exec)
    #:use-module (opencog bioscience)
    #:use-module (opencog ure)
    #:use-module (opencog pln)
    #:use-module (ice-9 threads)
    #:use-module (fibers conditions)
    #:use-module (ice-9 suspendable-ports)
    #:use-module (ice-9 textual-ports)
    #:use-module (srfi srfi-1)
    #:use-module (pln-bio bio-utils)
    #:use-module (pln-bio expr)
    #:use-module (pln-bio combo-preprocess)
    #:use-module (pln-bio preprocess))

(install-suspendable-ports!)

(define pred-var (Variable "$pred"))
(define gene-var (Variable "$gene"))
(define bp-var (Variable "$bp"))
(define num-rank 50)
(define batch-size 250)
(define sample-size 2237)
(define CT (Type "ConceptNode"))
(define GT (Type "GeneNode"))
(define PT (Type "PredicateNode"))
(define BT (TypeInh "ConceptNode"))

(define-public (run-concept-deduction-expr overexpr?)
    (setup-gene-expr overexpr?)
    ;; get patient atoms and run the deduction in batch
    (cog-logger-info "Generating SubsetLinks")
    ;;apply fc to get the relationship between go's and patients
    (let* ((patients (cog-get-atoms 'PatientNode))
            (batch-num 0)
            (batch-ls (split-lst patients batch-size))
            (batches (map (lambda (b) (set! batch-num (+ batch-num 1)) (cons batch-num b)) batch-ls))
            (prefix (if overexpr? "results/batches_overexpr_50/" "results/batches_underexpr_50/")))
        
        (if (not (file-exists? prefix)) (mkdir prefix))
        (n-par-for-each (/ (current-processor-count) 2)  (lambda (batch)
            (run-batch batch prefix overexpr?)) batches)
        (cog-logger-info "Done!")))

(define (run-batch batch prefix overexpr?)
    (let ((port (open-file (string-append prefix "batch_" (number->string (car batch)) ".scm") "w")))
        (for-each (lambda (patient)
            (if overexpr?
                (write-result-to-file port (generate-patient-concept-link-rule-overexpr patient))
                (write-result-to-file port (generate-patient-concept-link-rule-underexpr patient)))) (cdr batch))
        (close-port port)))

(define (generate-patient-concept-link-rule-overexpr patient-var)
    (cog-outgoing-set (cog-execute! (Bind 
        (VariableList
            (TypedVariable bp-var BT))
        (AndLink 
            (Present 
                patient-var
                bp-var)
            (SatisfactionLink 
                (TypedVariable (Variable "$gene") (Type "GeneNode"))
                (Present 
                    (EvaluationLink
                        (LazyExecutionOutputLink
                            (SchemaNode "make-overexpression-predicate-for-gene")
                            (VariableNode "$gene"))
                            patient-var)               
                    (Member (Variable "$gene") bp-var)
                    ;(Member (Variable "$gene") (Concept "top-ranked"))
                    )))
            
        (ExecutionOutputLink
            (GroundedSchemaNode "scm: generate-subset-tv-overexpr")
            (ListLink
                (Subset
                    (AndLink
                        bp-var
                        (Concept "profiled-genes"))              
                    (SatisfyingSetScope
                        (Variable "$G")
                        (EvaluationLink
                        (LazyExecutionOutputLink
                            (SchemaNode "make-overexpression-predicate-for-gene")
                            (VariableNode "$G"))
                            patient-var)))               
                patient-var
                bp-var))))))

(define (generate-patient-concept-link-rule-underexpr patient-var)
    (cog-outgoing-set (cog-execute! (Bind 
        (VariableList 
            (TypedVariable bp-var BT))
        (AndLink 
            (Present 
                patient-var
                bp-var)
            (SatisfactionLink 
                (TypedVariable (Variable "$gene") (Type "GeneNode"))
                (Present 
                    (EvaluationLink
                        (LazyExecutionOutputLink
                            (SchemaNode "make-underexpression-predicate-for-gene")
                            (VariableNode "$gene"))
                            patient-var)               
                    (Member (Variable "$gene") bp-var)
                    ;(Member (Variable "$gene") (Concept "top-ranked"))
                    )))
            
        (ExecutionOutputLink
            (GroundedSchemaNode "scm: generate-subset-tv-underexpr")
            (ListLink
                (Subset              
                    (AndLink 
                        bp-var
                        (Concept "profiled-genes"))
                    (SatisfyingSetScope
                        (Variable "$G")
                        (EvaluationLink
                        (LazyExecutionOutputLink
                            (SchemaNode "make-underexpression-predicate-for-gene")
                            (VariableNode "$G"))
                            patient-var)))               
                patient-var
                bp-var))))))

(define (create-lns-for-top-genes)
    (for-each (lambda (gene)
        (Member gene (Concept "top-ranked"))) (cog-get-atoms 'GeneNode)))

(define (inheritance->member)
    ;; Run FC to
    ;; 1. Translate Inheritance to MemberLink
    ;; 2. Infer closure of GO annotation

    (define (get-results-with-tvs result-lst) 
        (let ((all-mbrs (append (cog-outgoing-set result-lst)
            (get-member-links 'GeneNode 'ReactomeNode) (get-member-links 'GeneNode 'BiologicalProcessNode) 
            (get-member-links 'GeneNode 'MolecularFunctionNode) (get-member-links 'GeneNode 'CellularComponentNode) 
            (get-member-links 'GeneNode 'PharmGkbNode))))
            
            (filter all-nodes-non-null-mean? (map (lambda (x) (cog-set-tv! x (stv 1 1))) all-mbrs))))

        ;; Load PLN
        (pln-clear)
        (pln-load-from-file (get-full-path "rules/translation.scm"))
        (pln-load-from-file (get-full-path "rules/transitivity.scm"))
        (pln-add-rule 'present-inheritance-transitivity)
        (pln-add-rule 'present-mixed-member-inheritance-transitivity)
        (cog-logger-info "Running FC: inheritance->member")
        (write-atoms-to-file "concept-inhr-member.scm" (get-results-with-tvs (pln-fc source
            #:vardecl vardecl
            #:maximum-iterations mi
            #:complexity-penalty cp
            #:fc-full-rule-application fra))))

(define-public (generate-subset-tv-overexpr conclusion . premises)
    (if (= (length premises) 2)
        (let* ((patient (car premises))
               (go-term (cadr premises))
               (A-mbrs (get-profiled-go-mbrs go-term))
               (numt (numerator-overexpr A-mbrs patient))
               (dnt (denominator A-mbrs))
               (tvs (if (< 0 dnt) (exact->inexact (/ numt dnt)) 1))
               (tvc (count->confidence (length A-mbrs)))
               (TV (stv tvs tvc)))
               
            (if (< 0 tvc) (cog-merge-hi-conf-tv! conclusion TV)))))


(define-public (generate-subset-tv-underexpr conclusion . premises)
    (if (= (length premises) 2)
        (let* ((patient (car premises))
               (go-term (cadr premises))
               (A-mbrs (get-profiled-go-mbrs go-term))
               (numt (numerator-underexpr A-mbrs patient))
               (dnt (denominator A-mbrs))
               (tvs (if (< 0 dnt) (exact->inexact (/ numt dnt)) 1))
               (tvc (count->confidence (length A-mbrs)))
               (TV (stv tvs tvc)))
               
            (if (< 0 tvc) (cog-merge-hi-conf-tv! conclusion TV)))))

(define (get-profiled-go-mbrs A)
    (cog-outgoing-set (cog-execute! 
        (Get (TypedVariable (Variable "$g") (Type "GeneNode"))
        (Present 
            (Member 
                (Variable "$g")
                (Concept "profiled-genes"))           
            (Member 
                (Variable "$g")
                A))))))

;;Given the evaluation links of B and member links of A, calculate the fuzzy
;;intersection between patterns of A and B, expressed as

;; Sum_x min(pattern-of(X, A), pattern-of(X, B))

;; where pattern-of(X, A) is the strength of 
;; Member <TV>
;;  X
;;  A

;; whereas pattern-of(X, B) is the strength of 
;; Evaluation
;   (LazyExecutionOutput
;     (Schema "make-overexpression-predicate-for-gene")
;     X)
;   B))

(define-public (numerator-overexpr A-mbrs B)
    (define (fuzzy-intersect A-mbr)
        (let* ((execution-ln (cog-link 'LazyExecutionOutputLink (SchemaNode "make-overexpression-predicate-for-gene") A-mbr))
               (B-eval (cog-link 'EvaluationLink execution-ln B)))
            (if (null?  B-eval)
                0
                (min (cog-mean A-mbr) (cog-mean B-eval)))))

    (fold + 0 (map fuzzy-intersect A-mbrs)))


(define-public (numerator-underexpr A-mbrs B)
    (define (fuzzy-intersect A-mbr)
        (let* ((execution-ln (cog-link 'LazyExecutionOutputLink (SchemaNode "make-underexpression-predicate-for-gene") A-mbr))
               (B-eval (cog-link 'EvaluationLink execution-ln B)))
            (if (null?  B-eval)
                0
                (min (cog-mean A-mbr) (cog-mean B-eval)))))

    (fold + 0 (map fuzzy-intersect A-mbrs)))

;; Given evaluation links of B calculate the fuzzy sum of 
;;the patterns of B expressed as 
;; Sum_x pattern-of(X, B)
;; A-members

(define-public (denominator A-mbrs)
    (fold + 0 (map cog-mean A-mbrs)))

(define-public (find-go-genes-intersec go)
    (cog-outgoing-set (cog-execute! 
        (Get
            (VariableList 
                (TypedVariable (Variable "$g") (Type "GeneNode")))        
            (Present 
                (Member (Variable "$g") go)
                (Member (Variable "$g") (ConceptNode "profiled-genes")))))))

(define-public (setup-gene-expr overexpr?)
    ;; default opencog logger
    (cog-logger-set-stdout! #t)
    (cog-logger-set-filename! "logs/expr.log")
    ;;ure logger
    (ure-logger-set-timestamp! #f)
    (ure-logger-set-level! "debug")
    (ure-logger-set-filename! "logs/ure.log")

    (define filter-in (lambda (x)
        (or (go? x) (pathway? x)  (inheritance-GO? x) (pathway-inheritance? x)
            (gene-memberln-all? x))))

    ;;load go biological process
    (cog-logger-info "Loading knowlegde-base")
    (load-kbs (list "kbs/GO_2020-10-21.scm" "kbs/GO_annotation_2020-10-20.scm" 
    "kbs/reactome_2020-10-20.scm" "kbs/NCBI2Reactome_PE_Pathway.txt_2020-10-20.scm" "kbs/pharmagkb.scm" "kbs/concept-inhr-member.scm")
            #:filter-in filter-in)
    ;; Comment out after the first run and just load the output file
    ;(inheritance->member)
    (if overexpr?
        (begin 
            (cog-logger-info "Loading patient overexpression data")
            ;;load the atomese form of overexpr & underexpr
            (load-kbs (list "kbs/patient_gene_over_expr.scm"))

            ;;generate the quantiles for overexpr
            (cog-logger-info "Generating SchemaValueLists")
            (write-atoms-to-file "results/overexpr-dist_50.scm" (overexpression-dist))
            ;;get the evaluation links for overexpr
            (cog-logger-info "Generating EvaluationLinks")
            (write-atoms-to-file "results/overexpr-evals_50.scm" (get-overexpr-eval-ln))
            
            ;;Load moses models to get top ranked genes
            ;(cog-logger-info "Load moses models")
            ;(load-kbs (list "kbs/combo.scm"))
            ;;(create-lns-for-top-genes)
            )
        (begin 
            (cog-logger-info "Loading patient underexpression data")
            ;;load the atomese form of overexpr & underexpr
            (load-kbs (list "kbs/patient_gene_under_expr.scm"))

            ;;generate the quantiles for overexpr
            (cog-logger-info "Generating SchemaValueLists")
            (write-atoms-to-file "results/underexpr-dist_50.scm" (underexpression-dist))
            ;;get the evaluation links for overexpr
            (cog-logger-info "Generating EvaluationLinks")
            (write-atoms-to-file "results/underexpr-evals_50.scm" (get-underexpr-eval-ln))

            ;;Load moses models to get top ranked genes
            ;(cog-logger-info "Load moses models")
            ;(load-kbs (list "kbs/combo.scm"))
            ;(create-lns-for-top-genes)
        )))

(define-public (take-custom lst n)
    (if (< (length lst) n)
        (take lst (length lst))
        (take lst n)))

(define-public (drop-custom lst n)
    (if (< (length lst) n)
        (drop lst (length lst))
        (drop lst n)))
(define-public (split-lst lst n)
    (if (null? lst) '()
        (cons (take-custom lst n) (split-lst (drop-custom lst n) n))))