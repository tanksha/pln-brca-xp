(define-module (pln-bio expr)
    #:use-module (opencog)
    #:use-module (opencog exec)
    #:use-module (opencog bio)
)


(define personality_parm 10)
(define bin-size 3)


(define-public (overexpression-dist) 
    (cog-execute! (Bind
        (VariableList 
            (TypedVariable (Variable "$gene") (Type "GeneNode"))
            (TypedVariable (Variable "$patient") (Type "ConceptNode"))
            (TypedVariable (Variable "$num") (Type "NumberNode")))
        (Present 
            (ExecutionLink
                (ExecutionOutputLink
                    (SchemaNode "make-overexpression-schema-for-gene")
                    (Variable "$gene"))
                (Variable "$patient")
                (Variable "$num")))
                    
        (ExecutionOutputLink
            (GroundedSchema "scm: update-overexpr-svd")
            (List 
                (Variable "$gene")
                (Variable "$num")
                (SimpleTruthValue 1 1))))))



(define-public (update-overexpr-svd gene value overexpr)
        (cog-logger-debug "Updating svd for gene ~a:" gene)
        (let* ((overexpr? (stv->scm overexpr))
                (schema-record-ln (get-schema-val-record-ln gene #overexpr?)))
            (if (null? schema-record-ln)
                (if overexpr?
                    (SchemeValueRecordLink
                        (ExecutionOutputLink
                            (SchemaNode "make-overexpression-schema-for-gene")
                            gene)
                        (SchemaValueListLink
                            value
                            (NumberNode "1")))
                    (SchemeValueRecordLink
                        (ExecutionOutputLink
                            (SchemaNode "make-underexpression-schema-for-gene")
                            gene)
                        (SchemaValueListLink
                            value
                            (NumberNode "1"))))
                
                        
                (let* ((schema-lst (gdr (car schema-record-ln)))
                        (len (/ (length (cog-outgoing-set schema-lst)) 2)))
                        
                    (if (< len num-bins)
                        (update-schema-list-ln (car schema-record-ln) schema-lst gene value overexpr?)
                        (update-schema-record-ln (car schema-record-ln) schema-lst gene value overexpr?))))))


(define (get-schema-val-record-ln gene overexpr)
    (if overexpr
        (cog-outgoing-set (cog-execute! 
            (Bind 
                (SchemeValueRecordLink 
                    (ExecutionOutputLink
                        (SchemaNode "make-overexpression-schema-for-gene")
                        gene)
                    (Variable "$schema-lst"))
                (SchemeValueRecordLink 
                    (ExecutionOutputLink
                        (SchemaNode "make-overexpression-schema-for-gene")
                        gene)
                    (Variable "$schema-lst")))))
        (cog-outgoing-set (cog-execute! 
            (Bind 
                (SchemeValueRecordLink 
                    (ExecutionOutputLink
                        (SchemaNode "make-underexpression-schema-for-gene")
                        gene)
                    (Variable "$schema-lst"))
                (SchemeValueRecordLink 
                    (ExecutionOutputLink
                        (SchemaNode "make-underexpression-schema-for-gene")
                        gene)
                    (Variable "$schema-lst")))))))


(define (update-schema-record-ln record-ln schema-lst gene value overexpr?)
    (let* ((lst (cog-outgoin-set shema-lst))
           (closest (get-closet-value lst value))
           (v-cl (car closest))
           (n-cl (cadr closest))
           (index (cddr closest))
           (v-new (get-new-val value v-cl n-cl)))
           
        )
    )

(define (update-schema-lst-ln record-ln schema-lst gene value overexpr?)
    
    (let ((outgoing-nodes (append (cog-outgoing-set schema-lst) 
                (NumberNode value) (NumberNode "1")))) ;;what if a value occurs twice!
        ;;delete the schema record link then the schema list
        (cog-delete record-ln)
        (cog-delete schema-lst)
        (if overexpr?
            (SchemeValueRecordLink
                (SchemaNode "make-overexpression-schema-for-gene")
                (SchemaValueListLink outgoing-nodes))
            (SchemeValueRecordLink
                (SchemaNode "make-underexpression-schema-for-gene")
                (SchemaValueListLink outgoing-nodes)))))

(define (underexpression-dist) 
    (cog-execute! (Bind
        (VariableList 
            (TypedVariable (Variable "$gene") (Type "GeneNode"))
            (TypedVariable (Variable "$patient") (Type "ConceptNode"))
            (TypedVariable (Variable "$num") (Type "NumberNode")))
        (Present 
            (ExecutionLink
                (ExecutionLink
                    (SchemaNode "make-underexpression-schema-for-gene")
                    (Variable "$gene"))
                (Variable "$patient")
                (Variable "$num")))
                    
        (ExecutionOutputLink
            (GroundedSchema "scm: update-underexpr-svd")
            (List 
                (Variable "$gene")
                (Variable "$num"))))))

(define (quantile-helper svd i rem?)
    (let ((j (list-ref svd (- i 1)))
            (k (list-ref svd i)))
        (if rem? 
            (list j)
            (list (/ (+ k j) 2)))))

(define (get-quantile-borders lst q-size)
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
                    (append acc 
                        (if (= remainder 0) (quantile-helper lst i #f) '())
                        (list (list-ref lst (- size 1))))
                        
                    (loop (+ i step) (append acc (quantile-helper lst i #f)))))
              (else acc)))))

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

(define (get-new-val old-val vl-cl n-cl)
    ;;v_i_new = (v_i * n_i / (n_i + 1)) + (v / (n_i + 1))
    (+ (/ (* old-val n-cl) (+ n-cl 1)) (/ vl-cl (+ n-cl 1))))

(define (get-new-lst old-lst index val n-cl)
    (let lp ((i 0) (acc '()))
        (if (>= i (length old-lst) acc)
            (cond 
              ((= i index) (lp (+ i 1) (append acc (list val))))
              ((= i (+ index 1)) (lp (+ i 1) (append acc (list (+ n-cl 1)))))
              (else (lp (+ i 1) (append acc (list-ref old-lst i))))))))

(define (stv->scm tv)
  ;; do I need to include the confidence??
  (and (= 1 (cog-tv-mean tv)) (= 1 (cog-confidence tv))))