
(define-module (pln-bio preprocess)
    #:use-module (opencog)
    #:use-module (opencog exec)
    #:use-module (opencog logger)
    #:use-module (opencog ure)
    #:use-module (opencog pln)
    #:use-module (opencog bioscience)
    #:use-module (pln-bio bio-utils)
    #:use-module (pln-bio rules rule-utils)
    #:use-module (ice-9 threads)
    #:export (preprocess-int preprocess-ext)
)

;; Parameters
(define-public rs 0) ; Random seed
(define-public ss 0.3) ; Subsampled portion of the KBs
(define-public mi 12); Maximum number of iterations
(define-public cp 10); Complexity penalty
(define-public fra #t); Whether rules are fully applied

;; Defining this in a let scope causes segmentation fault
(define-public ConceptT (TypeInh "ConceptNode"))
(define-public GOT (TypeChoice go-types))
(define-public PathwayT (TypeChoice pathway-types))
(define-public GeneT (Type "GeneNode"))
(define-public X (Variable "$X"))
(define-public Y (Variable "$Y"))
(define-public vardecl (VariableSet
                  (TypedVariable X GOT)
                  (TypedVariable Y GOT)))
(define-public source (Inheritance X Y))

;;mutext to control file write access
(define mtx (make-mutex))

;; Parameters string
(define param-str (string-append
                   "-rs=" (number->string rs)
                   "-ss=" (number->string ss)
                   "-mi=" (number->string mi)
                   "-cp=" (number->string cp)
                   "-fra=" (bool->string fra)))


(define (inheritance->subset port)
    ;; Run FC to
    ;; 1. Translate Inheritance to SubSet
    ;; 2. Infer closure of GO annotation
    ;; Load PLN
    (cog-logger-info "Loading PLN rules")
    (pln-clear)
    (pln-load-from-path "pln-bio/rules/translation.scm")
    (pln-load-from-path "pln-bio/rules/transitivity.scm")
    (pln-add-rule 'present-inheritance-to-subset-translation)
    (pln-add-rule 'present-subset-transitivity)
    (pln-add-rule 'present-mixed-member-subset-transitivity)
    (cog-logger-info "PLN Rules loaded.")
    (cog-logger-info "Running FC: inheritance->subset")

    (n-par-for-each (current-processor-count) (lambda (x)
        (if (all-nodes-non-null-mean? x)
            (begin (cog-set-tv! x (stv 1 1)) 
                (if (lock-mutex mtx)
                    (begin (write x port) (unlock-mutex mtx))))))
        (cog-outgoing-set (pln-fc source
            #:vardecl vardecl
            #:maximum-iterations mi
            #:complexity-penalty cp
            #:fc-full-rule-application fra))))



(define (calculate-go/pathway-tvs all-atoms)
    ;; Calculate TVs of all GO categories or Pathways
    (define (concept-mean x sz) (exact->inexact (/ (get-cardinality x) sz)))
    (define (concept-tv x sz) (stv (concept-mean x sz) (count->confidence sz)))

    (let ((usize (length (get-genes))))
        (n-par-for-each (current-processor-count) (lambda (x) (cog-set-tv! x (concept-tv x usize))) all-atoms)))

(define (generate-subset port)
    (pln-clear)
    (pln-load-from-file "pln-bio/rules/subset-direct-introduction.scm")
    (pln-add-rule 'subset-direct-introduction)
    (define target (Subset X Y))
    (n-par-for-each (current-processor-count) (lambda (x)  
        (if (all-nodes-non-null-mean? x) 
            (if (lock-mutex mtx)
                (begin (write x port) (unlock-mutex mtx)))))                     
                            (cog-outgoing-set (pln-bc target #:vardecl vardecl
                                            #:maximum-iterations mi
                                            #:complexity-penalty cp))))

(define (generate-subset-negation port)
    (pln-clear)
    (pln-load-from-file "pln-bio/rules/subset-negation.scm")
    (pln-add-rule 'subset-condition-negation)
    (define target (Subset (Not X) Y))
    (n-par-for-each (current-processor-count) (lambda (x)
        (if (all-nodes-non-null-mean? x) 
           (if (lock-mutex mtx)
                (begin (write x port) (unlock-mutex mtx)))))  
                            (cog-outgoing-set (pln-bc target #:vardecl vardecl
                                            #:maximum-iterations mi
                                            #:complexity-penalty cp))))

(define (subset->attraction port)
   ;; Run backward chainer to produce attraction links. 
    ;; Add required PLN rules
    (pln-clear)
    (pln-load-from-file "pln-bio/rules/attraction-introduction.scm")
    (pln-add-rule 'subset-attraction-introduction)
    (define target (Attraction X Y))
    (n-par-for-each (current-processor-count) (lambda (x)
                (if (all-nodes-non-null-mean? x) 
                    (if (lock-mutex mtx)
                        (begin (write x port) (unlock-mutex mtx)))))                     
                            (cog-outgoing-set (pln-bc target #:vardecl vardecl
                                            #:maximum-iterations mi
                                            #:complexity-penalty cp))))


(define* (preprocess-int kbs #:key (filter-in #f))
   (let* ((scm-filename (string-append "results/sim/intensional/preprocess-kbs-asv2" param-str ".scm"))
         (port (open-file scm-filename "a")))

    ;;load kbs
    (cog-logger-info "Loading kbs")
    (if filter-in
        (load-kbs kbs #:subsmp ss #:filter-in filter-in)
        (load-kbs kbs #:subsmp ss))

    (cog-logger-info "Running FC: inheritance->subset")
    (inheritance->subset port)
    (cog-logger-info "Calculating GO Categories tvs")
    (calculate-go/pathway-tvs (get-go-categories))
    (cog-logger-info "Calculating Pathway tvs")
    ; (calculate-go/pathway-tvs (get-pathways) port)
    (cog-logger-info "Running BC: subset-introduction")
    (generate-subset port)
    (cog-logger-info "Running BC: subset-negation")
    (generate-subset-negation port)
    (cog-logger-info "Running BC: subset->attraction")
    (subset->attraction port)
    (cog-logger-info "Preprocessing done!")
    (close-port port)))

(define* (preprocess-ext kbs #:key (filter-in #f))
   (let* ((scm-filename (string-append "results/sim/extensional/preprocess-kbs-asv2" param-str ".scm"))
         (port (open-file scm-filename "a")))

    ;;load kbs
    (cog-logger-info "Loading kbs")
    (if filter-in
        (load-kbs kbs #:subsmp ss #:filter-in filter-in)
        (load-kbs kbs #:subsmp ss))

    (cog-logger-info "Running FC: inheritance->subset")
    (inheritance->subset port)
    (cog-logger-info "Calculating GO Categories tvs")
    (calculate-go/pathway-tvs (get-go-categories))
    (cog-logger-info "Calculating Pathway tvs")
    ; (calculate-go/pathway-tvs (get-pathways) port)
    (cog-logger-info "Preprocessing done!")
    (close-port port)))