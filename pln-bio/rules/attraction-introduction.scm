(define-module (pln-bio rules attraction-introduction)
   #:use-module (opencog) 
   #:use-module (opencog exec) 
   #:use-module (opencog ure) 
   #:use-module (opencog logger)
   #:use-module (srfi srfi-1)
   #:use-module (pln-bio rules rule-utils)
)

(define-public (gen-subset-attraction-introduction-rule A-types B-types)
  (let* ((A (Variable "$A"))
        (B (Variable "$B"))
        (A-type (TypeChoice A-types))
        (B-type (TypeChoice B-types)))
    (BindLink
      (VariableSet
        (TypedVariable A A-type)
        (TypedVariable B B-type))
      (Present
        (Subset A B)
        (Subset (Not A) B))
      (ExecutionOutputLink
        (GroundedSchemaNode "scm: attraction-introduction")
        (ListLink
          ;; Conclusion
          (Attraction A B)
          ;; Premises
          (Subset A B)
          (Subset (Not A) B))))))

;; Formula
(define (attraction-introduction conclusion . premises)
  (if (= (length premises) 2)
      (let* ((ATT conclusion)
             (SAB (car premises))
             (SNAB (cadr premises))
             (ATTs (max 0 (- (cog-mean SAB) (cog-mean SNAB))))
             (ATTc (min (cog-confidence SAB) (cog-confidence SNAB)))
             (ATTtv (stv ATTs ATTc)))
        (if (< 0 ATTc) (cog-merge-hi-conf-tv! ATT ATTtv)))))

; Name the rule

(define subset-attraction-introduction-rule
    (gen-subset-attraction-introduction-rule go-types go-types))

(define subset-attraction-introduction-rule-name
  (DefinedSchemaNode "subset-attraction-introduction-rule"))
(DefineLink subset-attraction-introduction-rule-name
  subset-attraction-introduction-rule)