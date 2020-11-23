(define-module (pln-bio expr)
    #:use-module (opencog)
    #:use-module (opencog exec)
    #:use-module (opencog bioscience)
    #:use-module (srfi srfi-1))


(define personality-parm 10)
(define bin-size 5)
(define num-samples 2237)
(define-public conf (exact->inexact (/ num-samples (+ num-samples personality-parm))))

(define-public (get-overexpr-eval-ln)
    (cog-execute! (Bind
        (VariableList
            (TypedVariable (Variable "$gene") (Type "GeneNode"))
            (TypedVariable (Variable "$patient") (Type "ConceptNode"))
            (TypedVariable (Variable "$num") (Type "NumberNode")))
        (Present 
            (ExecutionLink
                (LazyExecutionOutputLink
                    (SchemaNode "make-overexpression-schema-for-gene")
                    (Variable "$gene"))
                (Variable "$patient")
                (Variable "$num")))
        (ExecutionOutputLink 
            (GroundedSchemaNode "scm: calculate-expression-tv")
            (List 
                (Variable "$gene")
                (Variable "$patient")
                (Variable "$num")
                (Concept "1"))))))


(define-public (get-underexpr-eval-ln)
    (cog-execute! (Bind
        (VariableList
            (TypedVariable (Variable "$gene") (Type "GeneNode"))
            (TypedVariable (Variable "$patient") (Type "ConceptNode"))
            (TypedVariable (Variable "$num") (Type "NumberNode")))
        (Present 
            (ExecutionLink
                (LazyExecutionOutputLink
                    (SchemaNode "make-underexpression-schema-for-gene")
                    (Variable "$gene"))
                (Variable "$patient")
                (Variable "$num")))
        (ExecutionOutputLink 
            (GroundedSchemaNode "scm: calculate-expression-tv")
            (List 
                (Variable "$gene")
                (Variable "$patient")
                (Variable "$num")
                (ConceptNode "0"))))))

(define-public (calculate-expression-tv gene patient value overexpr-tv?)
    (let* ((overexpr? (= (num-node->num overexpr-tv?) 1))
          (quantiles (schema-list->nums (get-schema-lst gene overexpr?)))
          (tv (calculate-tv quantiles (num-node->num value))))
          
          (if overexpr?
            (Evaluation (stv tv conf)
                (LazyExecutionOutputLink
                    (SchemaNode "make-overexpression-predicate-for-gene")
                    gene)
                patient)
            (Evaluation (stv tv conf)
                (LazyExecutionOutputLink
                    (SchemaNode "make-underexpression-predicate-for-gene")
                    gene)
                patient))))

(define-public (overexpression-dist)
    (map (lambda (gene) 
        (create-schema-dist gene (get-overexpress-dist gene) #t)) (cog-get-atoms 'GeneNode)))

(define-public (underexpression-dist)
    (map (lambda (gene) 
        (create-schema-dist gene (get-underexpress-dist gene) #f)) (cog-get-atoms 'GeneNode)))

(define-public (create-schema-dist gene expr-lst overexpr?)
    (if (null? expr-lst) #f
        (let* ((sorted-vals (sort (get-values expr-lst) less))
            (quantiles (get-quantile-borders sorted-vals bin-size))
            (num-nodes (map (lambda (q) (Number q)) quantiles)))
        (if overexpr?
            (SchemaValueRecordLink
                (LazyExecutionOutputLink
                    (SchemaNode "make-overexpression-schema-for-gene")
                    gene)
                (SchemaValueListLink
                    num-nodes))
            (SchemaValueRecordLink
                (LazyExecutionOutputLink
                    (SchemaNode "make-underexpression-schema-for-gene")
                    gene)
                (SchemaValueListLink
                    num-nodes))))))


(define-public (get-overexpress-dist gene) 
    (cog-outgoing-set (cog-execute! (Bind
        (VariableList
            (TypedVariable (Variable "$patient") (Type "ConceptNode"))
            (TypedVariable (Variable "$num") (Type "NumberNode")))
        (Present 
            (ExecutionLink
                (LazyExecutionOutputLink
                    (SchemaNode "make-overexpression-schema-for-gene")
                    gene)
                (Variable "$patient")
                (Variable "$num"))
            (Member
                gene 
                (Concept "profiled-genes")))

        (ExecutionLink
                (LazyExecutionOutputLink
                    (SchemaNode "make-overexpression-schema-for-gene")
                    gene)
                (Variable "$patient")
                (Variable "$num"))))))

(define-public (get-underexpress-dist gene) 
    (cog-outgoing-set (cog-execute! (Bind
        (VariableList
            (TypedVariable (Variable "$patient") (Type "ConceptNode"))
            (TypedVariable (Variable "$num") (Type "NumberNode")))
        (Present 
            (ExecutionLink
                (LazyExecutionOutputLink
                    (SchemaNode "make-underexpression-schema-for-gene")
                    gene)
                (Variable "$patient")
                (Variable "$num"))
            (Member
                gene 
                (Concept "profiled-genes")))
                    
        (ExecutionLink
                (LazyExecutionOutputLink
                    (SchemaNode "make-underexpression-schema-for-gene")
                    gene)
                (Variable "$patient")
                (Variable "$num"))))))

(define-public (get-schema-lst gene overexpr?)
   (if overexpr? 
    (car (cog-outgoing-set (cog-execute! 
        (Get (SchemaValueRecordLink 
                (LazyExecutionOutputLink
                    (SchemaNode "make-overexpression-schema-for-gene")
                    gene)
                (Variable "$s"))))))
    (car (cog-outgoing-set (cog-execute! 
        (Get (SchemaValueRecordLink 
                (LazyExecutionOutputLink
                    (SchemaNode "make-underexpression-schema-for-gene")
                    gene)
                (Variable "$s"))))))))

(define-public (get-quantile-borders lst q-size)
    ;;Returns the the upper and lower bound  set of a quantile
    ;;eg. if svd =[10,20,30,40,50,60,70,80,90,100] and we need quartile,it will return
    ;;[10,25,45,65,100]

    (let* ((size (length lst))
            (step (euclidean-quotient size q-size))
            (remainder (euclidean-remainder size q-size)))
          
        (let loop ((i 0) (acc '()))
            (cond 
              ((= i 0) (loop (+ i step) (append acc (list (list-ref lst i)))))
              ((<= i size)
                (if (>= (+ i step) (- size 1))
                    (sort (delete-duplicates (append acc 
                        (if (= remainder 0) (quantile-helper lst i #f) '())
                        (list (list-ref lst (- size 1))))) less)
                        
                    (loop (+ i step) (append acc (quantile-helper lst i #f)))))
              (else (sort (delete-duplicates acc) less))))))

(define (quantile-helper svd i rem?)
    (let ((j (list-ref svd (- i 1)))
            (k (list-ref svd i)))
        (if rem? 
            (list j)
            (list (/ (+ k j) 2)))))

(define (calculate-quantile-weights quantiles)
    ;;given a list of quantiles it returns a list that contains
    ;;the quantiles and their weights
    (define size (length quantiles))
    (let lp ((i 0) (acc '()))
        (if (= i size) acc
            (lp (+ i 1) (append acc (list (list-ref quantiles i) (exact->inexact (/ i size))))))))

(define (get-closet-value lst val)

    (let ((closest-index 0)
          (diff 0) (size (length lst)))
          (let lp ((i 0)
                     (v-tmp (list-ref lst 0))
                     (diff-tmp (abs (- (list-ref lst 0) val))))
                (format #t "i - ~a, diff-tmp - ~a, diff - ~a\n" i diff-tmp diff)
                (if (= i 0) 
                    (begin 
                        (set! diff diff-tmp) 
                        (lp (+ i 2) v-tmp diff-tmp))
                    (if (< i size)
                        (if (< diff-tmp diff)
                            (begin 
                                (set! diff diff-tmp)
                                (set! closest-index i)
                                (lp (+ i 2) (list-ref lst i)
                                    (abs (- (list-ref lst i) val))))
                                
                            (lp (+ i 2) (list-ref lst i)
                                    (abs (- (list-ref lst i) val)))))))
                                
        (cons (list-ref lst closest-index) (cons (list-ref lst (+ closest-index 1)) closest-index))))

(define (less x y) (if (< x y) #t #f))

(define (get-values expr-lst)
        (map (lambda (expr)
            (string->number (cog-name (caddr (cog-outgoing-set expr))))) expr-lst))

(define (get-new-val old-val vl-cl n-cl)
    ;;v_i_new = (v_i * n_i / (n_i + 1)) + (v / (n_i + 1))
    (+ (/ (* old-val n-cl) (+ n-cl 1)) (/ vl-cl (+ n-cl 1))))

(define-public (calculate-tv lst v)
    (define size (length lst))
    (define qsize (- size 1))

    (define (calculate-tv-helper index)
        (let* ((k (list-ref lst (- index 1)))
              (j (list-ref lst index))
              (wk (exact->inexact (/ index qsize)))
              (wj (exact->inexact (/ (- index 1) qsize)))
              (diff (- j k)))
              
        (exact->inexact (/ (+ (* (- v k) wk) (* (- j v) wj)) diff))))

    (let lp ((i 0))
        (cond 
            ((and (= i 0) (<= v (list-ref lst i))) 0)
            ((and (> i 0) (<= v (list-ref lst i)))
                (calculate-tv-helper i))
            ((= i qsize) 1)
            (else (lp (+ i 1))))))

(define (get-closest-bin bins val)
    (define size (length bins))
    (let lp ((i 0))
        (if (>= i size) 0)
            (if (and (>= val (list-ref bins i)) (<= val (list-ref bins (+ i 2))))
                i
                (lp (+ i 2)))))

(define (get-new-lst old-lst index val n-cl)
    (let lp ((i 0) (acc '()))
        (if (>= i (length old-lst) acc)
            (cond 
              ((= i index) (lp (+ i 1) (append acc (list val))))
              ((= i (+ index 1)) (lp (+ i 1) (append acc (list (+ n-cl 1)))))
              (else (lp (+ i 1) (append acc (list-ref old-lst i))))))))

(define-public (schema-list->nums schema-ls)
    (map (lambda (node) (string->number (cog-name node)))  (cog-outgoing-set schema-ls)))

(define-public (num-node->num node) (string->number (cog-name node)))

(define (stv->scm tv)
  ;; do I need to include the confidence??
  (and (= 1 (cog-tv-mean tv)) (= 1 (cog-confidence tv))))