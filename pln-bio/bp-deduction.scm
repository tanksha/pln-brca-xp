(define-module (pln-bio bp-deduction)
    #:use-module (opencog)
    #:use-module (opencog exec)
    #:use-module (opencog bioscience)
    #:use-module (opencog ure)
    #:use-module (opencog pln)
    #:use-module (srfi srfi-1)
    #:use-module (pln-bio combo-preprocess))


(define pred-var (Variable "$pred"))
(define patient-var (Variable "$patient"))
(define gene-var (Variable "$gene"))
(define bp-var (Variable "$bp"))
(define num-rank 50)

(define CT (Type "ConceptNode"))
(define GT (Type "GeneNode"))
(define PT (Type "PredicateNode"))
(define BT (Type "BiologicalProcessNode"))

(define-public (generate-patient-bp-link-rule-overexpr)
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

(define-public (generate-patient-bp-link-rule-underexpr)
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

(define-public (create-lns-for-top-genes)
    (for-each (lambda (gene)
        (Member gene (Concept "top-ranked"))) (get-top-genes num-rank)))

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