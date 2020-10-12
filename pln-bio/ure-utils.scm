
(define-module (pln-bio ure-utils)
    #:use-module (opencog)
    #:use-module (opencog exec)
    #:use-module (opencog ure)
)

(define-public (get-direct-steps-to-target target)
"
  Return all inference steps directly inferring the given target, in
  the following format:

  (Set
    (List <rule-1> <source-1> <iteration-1>)
    ...
    (List <rule-n> <source-n> <iteration-n>))
"
  (let* ((pattern (Execution
                    (Variable "$rule")
                    (List
                      (Variable "$source")
                      (Variable "$iteration"))
                    target))
         (vardecl (VariableList
                    (TypedVariable (Variable "$rule") (Type 'DefinedSchemaNode))
                    (Variable "$source")
                    (TypedVariable (Variable "$iteration") (Type 'NumberNode))))
         (gl (Get vardecl pattern)))
    (cog-execute! gl)))

(define-public (get-direct-steps-from-source source)
"
  Return all inference steps directly inferred from the give source, in
  the following format:

  (Set
    (List <rule-1> <target-1> <iteration-1> )
    ...
    (List <rule-n> <target-n> <iteration-n>))
"
  (let* ((pattern (Execution
                    (Variable "$rule")
                    (List
                      source
                      (Variable "$iteration"))
                    (Variable "$target")))
         (vardecl (VariableList
                    (TypedVariable (Variable "$rule") (Type 'DefinedSchemaNode))
                    (Variable "$target")
                    (TypedVariable (Variable "$iteration") (Type 'NumberNode))))
         (gl (Get vardecl pattern)))
    (cog-execute! gl)))

(define-public (get-trails-to-target-rec target . inners)
"
  Return all inference trails leading to the given target, in the
  following format:

  (Set
    (List
      (List <rule-11> <inter-11> <iteration-11>)
      ...
      (List <rule-1m> <inter-1m> <iteration-1m>))
    ...
    (List
      (List <rule-n1> <inter-n1> <iteration-n1>)
      ...
      (List <rule-nm> <inter-nm> <iteration-nm>)))
"
  (let* ((get-inner (lambda (s) (gdr s))) ; Get the inner target of a step
         (direct-steps (get-direct-steps-to-target target))
         ;; Remove cycles
         (inners? (lambda (s) (member (get-inner s) inners)))
         (not-inners? (lambda (s) (not (inners? s))))
         (direct-steps-no-cycles (filter not-inners? (cog-outgoing-set direct-steps)))
         ;; Given a direct inference step, find the trails going to
         ;; that inference step, and append the inference step to them
         (get-trails (lambda (s)
                       (let* ((inrs (if (inners? s) inners (cons (get-inner s) inners))))
                         (cog-outgoing-set (apply get-trails-to-target-rec (cons (get-inner s) inrs))))))
         (append-step-to-trail (lambda (t s)
                                 (List (cog-outgoing-set t) s)))
         (append-step-to-trails (lambda (ts s)
                                  (if (null? ts)
                                      (List s)
                                      (map (lambda (t) (append-step-to-trail t s)) ts))))
         (get-trails-with-direct-step (lambda (s)
                                        (let* ((ts (get-trails s)))
                                          (append-step-to-trails ts s)))))
    (Set (map get-trails-with-direct-step direct-steps-no-cycles))))

(define-public (get-trails-to-target target)
  (get-trails-to-target-rec target target))

(define-public (get-trails-from-source-rec source . inners)
"
  Return all inference trails coming from the given source, in the
  following format:

  (Set
    (List
      (List <rule-11> <inter-11> <iteration-11>)
      ...
      (List <rule-1m> <inter-1m> <iteration-1m>))
    ...
    (List
      (List <rule-n1> <inter-n1> <iteration-n1>)
      ...
      (List <rule-nm> <inter-nm> <iteration-nm>)))
"
  (let* ((get-inner (lambda (s) (gdr s))) ; Get the inner target of a step
         (direct-steps (get-direct-steps-from-source source))
         ;; Remove cycles
         (inners? (lambda (s) (member (get-inner s) inners)))
         (not-inners? (lambda (s) (not (inners? s))))
         (direct-steps-no-cycles (filter not-inners? (cog-outgoing-set direct-steps)))
         ;; Given a direct inference step, find the trails going to
         ;; that inference step, and append the inference step to them
         (get-trails (lambda (s)
                       (let* ((inrs (if (inners? s) inners (cons (get-inner s) inners))))
                         (cog-outgoing-set (apply get-trails-from-source-rec (cons (get-inner s) inrs))))))
         (prepend-step-to-trail (lambda (t s) (List s (cog-outgoing-set t))))
         (prepend-step-to-trails (lambda (ts s)
                                   (if (null? ts)
                                       (List s)
                                       (map (lambda (t) (prepend-step-to-trail t s)) ts))))
         (get-trails-with-direct-step (lambda (s)
                                        (let* ((ts (get-trails s)))
                                          (prepend-step-to-trails ts s)))))
    (Set (map get-trails-with-direct-step direct-steps-no-cycles))))

(define-public (get-trails-from-source source)
  (get-trails-from-source-rec source source))