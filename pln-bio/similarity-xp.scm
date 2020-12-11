(define-module (pln-bio similarity-xp)
    #:use-module (opencog)
    #:use-module (opencog exec)
    #:use-module (opencog randgen)
    #:use-module (opencog logger)
    #:use-module (opencog ure)
    #:use-module (opencog pln)
    #:use-module (opencog bioscience)
    #:use-module (pln-bio bio-utils)
    #:use-module (pln-bio preprocess)
)

(define X (Variable "$X"))
(define Y (Variable "$Y"))
(define ConceptT (TypeInh "ConceptNode"))
(define int-target (IntensionalSimilarity X Y))
(define ext-target (ExtensionalSimilarity X Y))
(define vardecl (VariableSet
                  (TypedVariable X ConceptT)
                  (TypedVariable Y ConceptT)))

(define-public (go-pathway-intensional-similarity kbs)
   (define opencog-log-filename "logs/extensional-reasoning-test.log")
   (define ure-log-filename "logs/extensional-reasoning-test-ure.log")

    (cog-logger-set-timestamp! #f)
    (cog-logger-set-stdout! #t)
    (cog-logger-set-sync! #t)
    (cog-logger-set-level! "info")
    (cog-logger-set-filename! opencog-log-filename)
    (ure-logger-set-timestamp! #f)
    (ure-logger-set-sync! #t)
    (ure-logger-set-level! "debug")
    (ure-logger-set-filename! ure-log-filename)

    (let* ((rs 0) (mi 100) (cp 1)
           (param-str (string-append
                   "-rs=" (number->string rs)
                   "-mi=" (number->string mi)
                   "-cp=" (number->string cp)))
            (output-file (string-append "results/sim/go-intensional-similarity" param-str ".scm"))
            (filter-in (lambda (x)
                            (or (go? x)  (inheritance-GO_term? x)
                                (gene-memberln? x)))))

        ;;Preprocessing step
        (preprocess kbs #:filter-in filter-in)
        ;; Load PLN
        (cog-logger-info "Running BC: Attraction->IntensionalSimilarity")
        (pln-clear)
        (pln-load-from-path "opencog/pln/rules/intensional/intensional-similarity-direct-introduction.scm")
        (pln-add-rule 'intensional-similarity-direct-introduction)

        (write-atoms-to-file output-file (cog-outgoing-set (pln-bc int-target 
            #:vardecl vardecl
            #:maximum-iterations mi #:complexity-penalty cp)))
        (cog-logger-info "Done!")))

(define-public (go-pathway-extensional-similarity kbs)
   (define opencog-log-filename "logs/extensional-reasoning-test.log")
   (define ure-log-filename "logs/extensional-reasoning-test-ure.log")

    (cog-logger-set-timestamp! #f)
    (cog-logger-set-stdout! #t)
    (cog-logger-set-sync! #t)
    (cog-logger-set-level! "info")
    (cog-logger-set-filename! opencog-log-filename)
    (ure-logger-set-timestamp! #f)
    (ure-logger-set-sync! #t)
    (ure-logger-set-level! "debug")
    (ure-logger-set-filename! ure-log-filename)

    (let* ((rs 0) (mi 100) (cp 1)
           (param-str (string-append
                   "-rs=" (number->string rs)
                   "-mi=" (number->string mi)
                   "-cp=" (number->string cp)))
            (output-file (string-append "results/sim/go-extensional-similarity" param-str ".scm"))
            (filter-in (lambda (x)
                            (or (go? x) (inheritance-GO_term? x)
                                (gene-memberln? x)))))

        ;;Preprocessing step
        (preprocess kbs #:filter-in filter-in)
        ;; Load PLN
        (cog-logger-info "Running BC: Attraction->ExtensionalSimilarity")
        (pln-clear)
        (pln-load-from-path "opencog/pln/rules/extensional/extensional-similarity-direct-introduction.scm")
        (pln-add-rule 'extensional-similarity-direct-introduction)

        (write-atoms-to-file output-file (cog-outgoing-set (pln-bc ext-target 
            #:vardecl vardecl
            #:maximum-iterations mi #:complexity-penalty cp)))
        (cog-logger-info "Done!")))