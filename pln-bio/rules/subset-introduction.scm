(define-module (pln-bio rules subset-introduction)
   #:use-module (opencog) 
   #:use-module (opencog exec) 
   #:use-module (opencog ure) 
   #:use-module (opencog logger)
   #:use-module (srfi srfi-1)
   #:use-module (pln-bio rules rule-utils)
   #:use-module (pln-bio rules extensional-utils)
)

;; Rule generator for subset introduction with specific type
(define-public (gen-subset-direct-introduction A-types B-types)
    (let* ((A (Variable "$A"))
           (B (Variable "$B"))
           (A-type (TypeChoice A-types (Type "AndLink")))
           (B-type (TypeChoice B-types (Type "AndLink"))))
        
        (Bind 
            (VariableSet 
                (TypedVariable A A-type)
                (TypedVariable B B-type))
            (Present
                A
                B)
            (ExecutionOutput
                (GroundedSchema "scm: subset-direct-introduction")
                (List
                ;; Conclusion
                (Subset A B)
                ;; Premises
                A
                B)))))

;; Given a list of members of A and B calculate the TV of
;;
;; SubsetLink
;;   A
;;   B
(define (subset-evidence->tv A-mbrs B-mbrs)
  ;; (cog-logger-debug "(subset-evidence->tv A-mbrs=~a B-mbrs=~a)" A-mbrs B-mbrs)
  (let* ;; TODO consider TVs of the members
       ((A-size (length A-mbrs))
        (AB-mbrs (lset-intersection equal? A-mbrs B-mbrs))
        (AB-size (length AB-mbrs))
        (strength (if (< 0 A-size)
                      (exact->inexact (/ AB-size A-size))
                      1))
        (confidence (if (< 0 A-size)
                        (count->confidence A-size)
                        0)))
    (stv strength confidence)))

;; Formula
(define (subset-direct-introduction conclusion . premises)
  ;; (cog-logger-debug "(subset-direct-introduction conclusion=~a . premises=~a)" conclusion premises)
  (if (= (length premises) 2)
      (let* ((Ss conclusion)
             (A (car premises))
             (B (cadr premises))
             ;; Fetch all members of A and B
             (A-mbrs (get-members-of A))
             (B-mbrs (get-members-of B))
             ;; Calculate the TV based on the members of A and B
             (tv (subset-evidence->tv A-mbrs B-mbrs)))
        (if (and (< 0 (cog-tv-mean tv)) (< 0 (cog-tv-confidence tv)))
            (cog-merge-hi-conf-tv! Ss tv)))))

(define subset-direct-introduction-rule
    (gen-subset-direct-introduction go-types go-types))

(define subset-direct-introduction-rule-name
  (DefinedSchemaNode "subset-direct-introduction-rule"))
(DefineLink subset-direct-introduction-rule-name
  subset-direct-introduction-rule)