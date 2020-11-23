(define-module (pln-bio bp-deduction)
    #:use-module (opencog)
    #:use-module (opencog logger)
    #:use-module (opencog exec)
    #:use-module (opencog bioscience)
    #:use-module (opencog ure)
    #:use-module (opencog pln)
    #:use-module (srfi srfi-1)
    #:use-module (pln-bio bio-utils)
    #:use-module (pln-bio expr)
    #:use-module (pln-bio combo-preprocess)
    #:use-module (pln-bio preprocess))


(define pred-var (Variable "$pred"))
(define patient-var (Variable "$patient"))
(define gene-var (Variable "$gene"))
(define bp-var (Variable "$bp"))
(define num-rank 50)

(define CT (Type "ConceptNode"))
(define GT (Type "GeneNode"))
(define PT (Type "PredicateNode"))
(define BT (Type "BiologicalProcessNode"))


(define (generate-patient-bp-link-rule-overexpr)
    (cog-execute! (Bind 
        (VariableList 
            (TypedVariable patient-var CT)
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
                    (Member (Variable "$gene") (Concept "top-ranked")))))
            
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
                bp-var)))))

(define (generate-patient-bp-link-rule-underexpr)
    (cog-execute! (Bind 
        (VariableList 
            (TypedVariable patient-var CT)
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
                    (Member (Variable "$gene") (Concept "top-ranked")))))
            
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
                bp-var)))))

(define (create-lns-for-top-genes)
    (for-each (lambda (gene)
        (Member gene (Concept "top-ranked"))) (get-top-genes num-rank)))

(define (inheritance->member)
    ;; Run FC to
    ;; 1. Translate Inheritance to MemberLink
    ;; 2. Infer closure of GO annotation

    (define (get-results-with-tvs result-lst) 
        (let ((all-mbrs (append (cog-outgoing-set result-lst)
                    (get-member-links 'GeneNode 'BiologicalProcessNode))))
            
            (filter all-nodes-non-null-mean? (map (lambda (x) (cog-set-tv! x (stv 1 1))) all-mbrs))))

        ;; Load PLN
        (pln-clear)
        (pln-load-from-file (get-full-path "rules/translation.scm"))
        (pln-load-from-file (get-full-path "rules/transitivity.scm"))
        (pln-add-rule 'present-inheritance-transitivity)
        (pln-add-rule 'present-mixed-member-inheritance-transitivity)
        (cog-logger-info "Running FC: inheritance->member")
        (write-atoms-to-file "go-inhr-member.scm" (get-results-with-tvs (pln-fc source
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

(define-public (run-expr-deduction overexpr?)
    ;; default opencog logger
    (cog-logger-set-stdout! #t)
    (cog-logger-set-filename! "logs/expr.log")
    ;;ure logger
    (ure-logger-set-timestamp! #f)
    (ure-logger-set-level! "debug")
    (ure-logger-set-filename! "logs/ure.log")

    (define filter-in (lambda (x)
        (or (go_bp? x)  (inheritance-GO_bp? x)
            (gene-memberln? x))))

    ;;load go biological process
    (cog-logger-info "Loading GO terms")
    (load-kbs (list "kbs/GO.scm" "kbs/GO_annotation.scm")
            #:filter-in filter-in)
    ;; Comment out after the first run and just load the output file
    (inheritance->member)
    (if overexpr?
        (begin 
            (cog-logger-info "Loading patient overexpression data")
            ;;load the atomese form of overexpr & underexpr
            (load-kbs (list "kbs/patient_gene_over_expr_50genes.scm"))

            ;;generate the quantiles for overexpr
            (cog-logger-info "Generating SchemaValueLists")
            (write-atoms-to-file "results/overexpr-dist.scm" (overexpression-dist))
            ;;get the evaluation links for overexpr
            (cog-logger-info "Generating EvaluationLinks")
            (write-atoms-to-file "results/overexpr-evals.scm" (cog-outgoing-set (get-overexpr-eval-ln)))
            
            ;;Load moses models to get top ranked genes
            (cog-logger-info "Load moses models")
            (load-kbs (list "kbs/combo.scm"))
            (create-lns-for-top-genes)

            (cog-logger-info "Generating SubsetLinks")
            ;;apply fc to get the relationship between go's and patients
            (write-atoms-to-file "results/subset-bp-patient-overexpr_50genes.scm" (cog-outgoing-set (generate-patient-bp-link-rule-overexpr)))
            (cog-logger-info "Done!")
        )
        (begin 
            (cog-logger-info "Loading patient underexpression data")
            ;;load the atomese form of overexpr & underexpr
            (load-kbs (list "kbs/patient_gene_under_expr_50genes.scm"))

            ;;generate the quantiles for overexpr
            (cog-logger-info "Generating SchemaValueLists")
            (write-atoms-to-file "results/underexpr-dist.scm" (underexpression-dist))
            ;;get the evaluation links for overexpr
            (cog-logger-info "Generating EvaluationLinks")
            (write-atoms-to-file "results/underexpr-evals.scm" (cog-outgoing-set (get-underexpr-eval-ln)))

            ;;Load moses models to get top ranked genes
            (cog-logger-info "Load moses models")
            (load-kbs (list "kbs/combo.scm"))
            (create-lns-for-top-genes)

            (cog-logger-info "Generating SubsetLinks")
            ;;apply fc to get the relationship between go's and patients
            (pln-clear)
            (write-atoms-to-file "results/subset-bp-patient-underexpr_50genes.scm" (cog-outgoing-set (generate-patient-bp-link-rule-underexpr)))
            (cog-logger-info "Done!"))))