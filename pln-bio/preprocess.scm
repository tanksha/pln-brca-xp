
(define-module (pln-bio preprocess)
    #:use-module (opencog)
    #:use-module (opencog exec)
    #:use-module (opencog logger)
    #:use-module (opencog ure)
    #:use-module (opencog pln)
    #:use-module (opencog bioscience)
    #:use-module (pln-bio bio-utils)
    #:use-module (pln-bio rule-utils)
    #:export (preprocess)
)

;; Parameters
(define-public rs 0) ; Random seed
(define-public ss 1) ; Subsampled portion of the KBs
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
;; Parameters string
(define param-str (string-append
                   "-rs=" (number->string rs)
                   "-ss=" (number->string ss)
                   "-mi=" (number->string mi)
                   "-cp=" (number->string cp)
                   "-fra=" (bool->string fra)))


(define (inheritance->subset)
    ;; Run FC to
    ;; 1. Translate Inheritance to SubSet
    ;; 2. Infer closure of GO annotation

    (define (get-results-with-tvs result-lst) 
        (let ((all-mbrs (append (cog-outgoing-set result-lst)
                    (get-member-links 'GeneNode 'ReactomeNode)
                    (get-member-links 'GeneNode 'SmpNode)
                    (get-member-links 'GeneNode 'PharmGkbNode)
                    (get-member-links 'GeneNode 'BiologicalProcessNode))))
            
            (filter all-nodes-non-null-mean? (map (lambda (x) (cog-set-tv! x (stv 1 1))) all-mbrs))))

        ;; Load PLN
        (cog-logger-info "Loading PLN rules")
        (pln-clear)
        (pln-load-from-file (get-full-path "rules/translation.scm"))
        (pln-load-from-file (get-full-path "rules/transitivity.scm"))
        (pln-add-rule 'present-inheritance-to-subset-translation)
        (pln-add-rule 'present-subset-transitivity)
        (pln-add-rule 'present-mixed-member-subset-transitivity)
        (cog-logger-info "PLN Rules loaded.")
        (cog-logger-info "Running FC: inheritance->subset")
        (get-results-with-tvs (pln-fc source
            #:vardecl vardecl
            #:maximum-iterations mi
            #:complexity-penalty cp
            #:fc-full-rule-application fra)))



(define (calculate-go/pathway-tvs all-atoms)
    ;; Calculate TVs of all GO categories or Pathways
    (define (concept-mean x sz) (exact->inexact (/ (get-cardinality x) sz)))
    (define (concept-tv x sz) (stv (concept-mean x sz) (count->confidence sz)))

    (let ((usize (length (get-genes))))
        (filter all-nodes-non-null-mean? (map (lambda (x) (cog-set-tv! x (concept-tv x usize))) all-atoms))))

(define (get-inverse-subsets atoms)
    ; Infer inverse subset links for pathways or go categories, based on inversion
    (filter (lambda (x) (and (gt-zero-mean-and-confidence? x) (all-nodes-non-null-mean? x))) (map true-subset-inverse atoms)))



(define (generate-subset)
    (pln-clear)
    (pln-load-from-file (get-full-path "rules/subset-direct-introduction.scm"))
    (pln-add-rule 'subset-direct-introduction)
    (define target (Subset X Y))
    (filter all-nodes-non-null-mean? (cog-outgoing-set (pln-bc target #:vardecl vardecl
                                            #:maximum-iterations mi
                                            #:complexity-penalty cp))))

(define (generate-subset-negation)
    (pln-clear)
    (pln-load-from-file (get-full-path "rules/subset-negation.scm"))
    (pln-add-rule 'subset-condition-negation)
    (define target (Subset (Not X) Y))
    (filter all-nodes-non-null-mean? (cog-outgoing-set (pln-bc target #:vardecl vardecl
                                            #:maximum-iterations mi
                                            #:complexity-penalty cp)))
)

(define (subset->attraction)
   ;; Run backward chainer to produce attraction links. 
    ;; Add required PLN rules
    (pln-clear)
    (pln-load-from-file (get-full-path "rules/attraction-introduction.scm"))
    (pln-add-rule 'subset-attraction-introduction)
    (define target (Attraction X Y))
    (filter all-nodes-non-null-mean? (cog-outgoing-set (pln-bc target #:vardecl vardecl
                                            #:maximum-iterations mi
                                            #:complexity-penalty cp))))


(define* (preprocess kbs #:key (filter-in #f))
   (let ((scm-filename (string-append "results/sim/preprocess-kbs-asv2" param-str ".scm"))
         (port (open-file scm-file "a")))

    ;;load kbs
    (cog-logger-info "Loading kbs")
    (if filter-in
        (load-kbs kbs #:subsmp ss #:filter-in filter-in)
        (load-kbs kbs #:subsmp ss))

    (write-result-to-file port (inheritance->subset))
    (cog-logger-info "Calculating GO Categories tvs")
    (write-result-to-file port (calculate-go/pathway-tvs (get-go-categories)))
    (cog-logger-info "Calculating Pathway tvs")
    ; (write-result-to-file port (calculate-go/pathway-tvs (get-pathways)))
    (cog-logger-info "Running BC: subset-introduction")
    (write-result-to-file port (generate-subset))
    (cog-logger-info "Running BC: subset-negation")
    (write-result-to-file port (generate-subset-negation))
    ; (cog-logger-info "Getting inverse Pathway Subsets")
    ; (write-result-to-file port (get-inverse-subsets (get-pathway-subsets)))
    (cog-logger-info "Running BC: subset->attraction")
    (write-result-to-file port (subset->attraction))
    (cog-logger-info "Preprocessing done!")
    (close-port port)))