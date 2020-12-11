(use-modules (pln-bio rule-utils))

(define (gen-subset-condition-negation-rule A-types B-types)
  (let* ((A (Variable "$A"))
        (B (Variable "$B"))
        (A-type (TypeChoice A-types))
        (B-type (TypeChoice B-types)))
    (Bind
      (VariableSet
        (TypedVariable A A-type)
        (TypedVariable B B-type))
      (Present
        (Subset A B))
      (ExecutionOutput
        (GroundedSchema "scm: subset-condition-negation")
        (List
          ;; Conclusion
          (Subset (Not A) B)
          ;; Premises
          (Subset A B)
          A
          B)))))

;; Formula
(define (subset-condition-negation conclusion . premises)
  (if (= (length premises) 3)
      (let* ((NS conclusion)
             (S (car premises))
             (A (cadr premises))
             (B (caddr premises))
             (Ss (cog-mean S))
             (As (cog-mean A))
             (Ac (cog-confidence A))
             (Bs (cog-mean B))
             (NSs (if (< As 1)
                      (/ (- Bs (* Ss As)) (- 1 As))
                      1))
             (NSc (if (< As 1) Ac 0))
             (NStv (stv NSs NSc)))
        (cog-merge-hi-conf-tv! NS NStv))))

(define subset-condition-negation-rule
    (gen-subset-condition-negation-rule go-types go-types))

;; Name
(define subset-condition-negation-rule-name
  (DefinedSchemaNode "subset-condition-negation-rule"))
(DefineLink subset-condition-negation-rule-name subset-condition-negation-rule)