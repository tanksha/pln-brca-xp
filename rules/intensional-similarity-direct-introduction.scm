(use-modules (pln-bio rule-utils))

(define (gen-intensional-similarity-direct-introduction-rule A-types B-types)
  (let* ((A (Variable "$A"))
           (B (Variable "$B"))
           (A-type (TypeChoice A-types))
           (B-type (TypeChoice B-types)))
    (Bind
      (VariableSet
        (TypedVariable A A-type)
        (TypedVariable B B-type))
      (And
        (Present
          A
          B)
        ;; There exists X such that
        ;;
        ;; (Attraction A X)
        ;; (Attraction B X)
        ;;
        ;; are present in the atomspace
        (Satisfaction
          (TypedVariable X CT)
          (Present
            (Attraction A X)
            (Attraction B X)))
        ;; A and B are different
        (Not (Equal A B)))
      (ExecutionOutput
        (GroundedSchema "scm: intensional-similarity-direct-introduction")
        (List
          ;; Conclusion
          (IntensionalSimilarity A B)
          ;; Premises (wrapped in Set because commutative)
          (Set  
            A
            B))))))

;; Formula
(define (intensional-similarity-direct-introduction conclusion . premises)
  ;; Given a concept return all attraction link
  ;;
  ;; Attraction <TV>
  ;;   A
  ;;   X
  (define (get-attractions A)
    (let* ((at-links (cog-filter 'AttractionLink (cog-incoming-set A)))
           (A-at? (lambda (x) (equal? A (gar x)))))
      (filter A-at? at-links)))

  ;; The pattern strength is the product of the mean and the
  ;; confidence of the TV on the attraction link
  ;;
  ;; Attraction <TV>
  ;;   A
  ;;   X
  (define (get-pattern-strength A pat)
    (let* ((A-at (cog-link 'AttractionLink A pat)))      
      (if (null? A-at) 0 (* (cog-mean A-at) (cog-confidence A-at)))))

  ;; Given the attraction links of A and B calculate the fuzzy
  ;; intersection between the patterns of A and B, expressed as
  ;;
  ;; Sum_x min(pattern-of(X,A), pattern-of(X,B))
  ;;
  ;; where pattern-of(X,A) is the strength of the TV of
  ;;
  ;; Attraction <TV>
  ;;   A
  ;;   X
  (define (numerator A B-ats)
    (define (fuzzy-intersect B-at)
      (let* ((pat (gdr B-at)))
        (min (get-pattern-strength A pat) (cog-mean B-at))))
    (fold + 0 (map fuzzy-intersect B-ats)))

  ;; Given the attraction links of A and B calculate the fuzzy sum of
  ;; the union of patterns of A and B expressed as
  ;;
  ;; Sum_x max(pattern-of(X,A), pattern-of(X,B))
  (define (denominator A B pats)
    (define (fuzzy-union pat)
      (let ((A-pat-strength (get-pattern-strength A pat))
            (B-pat-strength (get-pattern-strength B pat)))
        (max A-pat-strength B-pat-strength)))
    (fold + 0 (map fuzzy-union pats)))

  ;; (cog-logger-debug "(intensional-similarity-direct-introduction conclusion=~a . premises=~a)" conclusion premises)
  (if (= (length premises) 1)
      (let* ((IntInh conclusion)
             (A (gar (car premises)))
             (B (gdr (car premises)))
             ;; Fetch all pattern attraction links and patterns
             (A-ats (get-attractions A))
             (B-ats (get-attractions B))
             (A-pats (map gdr A-ats))
             (B-pats (map gdr B-ats))
             (pats (lset-union equal? A-pats B-pats))
             ;; Calculate denominator, then
             (dnt (denominator A B pats))
             (TVs (if (< 0 dnt) (/ (numerator A B-ats) dnt) 1))
             (TVc (count->confidence (length B-ats)))
             (TV (stv TVs TVc)))
        (if (< 0 TVc) (cog-merge-hi-conf-tv! IntInh TV)))))

; Name the rule
(define intensional-similarity-direct-introduction-rule (gen-intensional-similarity-direct-introduction-rule go-types go-types))

(define intensional-similarity-direct-introduction-rule-name
  (DefinedSchemaNode "intensional-similarity-direct-introduction-rule"))
(DefineLink intensional-similarity-direct-introduction-rule-name
  intensional-similarity-direct-introduction-rule)