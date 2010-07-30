
(ns chemiclj.smiles2
  (:use [chemiclj.core])
  (:require [shortcut.graph :as graph]
            [edu.arizona.fnparse [hound :as h] [core :as c]]
            [clojure.contrib [except :as except]])
  (:refer-clojure :exclude #{read-string}))

(defrecord BracketAtomContext
  [isotope symbol chirality hydrogen-count charge class])

(defn str* [objects]
  (apply str objects))

(h/defrule <left-bracket> (h/lit \[))
(h/defrule <right-bracket> (h/lit \]))

(h/defrule <aliphatic-element-symbol>
  (h/hook (comp str* concat)
          (h/cat h/<uppercase-ascii-letter>
                 (h/opt h/<lowercase-ascii-letter>)
                 (h/opt h/<lowercase-ascii-letter>))))

(h/defrule <aromatic-element-symbol>
  (h/+ (h/hook (comp str* concat) (h/cat (h/lit \s) (h/lit \e)))
       (h/hook (comp str* concat) (h/cat (h/lit \a) (h/lit \s)))
       (h/hook str (h/lit \c))
       (h/hook str (h/lit \n))
       (h/hook str (h/lit \o))
       (h/hook str (h/lit \p))
       (h/hook str (h/lit \s))))

(h/defrule <element-symbol>
  (h/label "an element symbol"
           (h/+ <aliphatic-element-symbol> <aromatic-element-symbol>)))

(h/defmaker radix-natural-number [core]
  (h/hooked-rep #(+ (* core %1) %2) 0 (h/radix-digit core)))

(h/defrule <decimal-natural-number>
  (radix-natural-number 10))

(h/defrule <isotope>
  (h/label "an isotope" <decimal-natural-number>))

(h/defrule <bracket-mods>
  (h/hook (fn [x] (when (seq x) (apply conj x)))
          (h/rep*
           (h/+
            (h/hook (fn [[_ hydrogen-count]]
                      {:hydrogens (or hydrogen-count 1)})
                    (h/cat (h/lit \H) (h/opt <decimal-natural-number>)))
            (h/hook (fn [l] {:chirality (str* (concat l))})
                    (h/+ (h/cat (h/lit \@)) (h/lit \@)
                         (h/lit \@)))
            (h/hook (fn [[sgn num]]
                      {:charge ((if (= sgn \+) + -) num)})
                    (h/+ (h/cat (h/lit \+) (h/opt <decimal-natural-number>))
                         (h/cat (h/lit \-) (h/opt <decimal-natural-number>))))))))

(h/defrule <bracket-expr>
  (h/hook (fn [[isotope symbol {:keys #{hydrogens chirality charge}} class]]
            (chemiclj.smiles2.BracketAtomContext.
             isotope symbol chirality hydrogens charge class))
          (h/label "a bracket expression"
                   (h/prefix <left-bracket>
                             (h/suffix
                              (h/cat
                               (h/opt <isotope>)
                               <element-symbol>
                               <bracket-mods>
                               (h/opt
                                (h/prefix (h/lit \:)
                                          <decimal-natural-number>)))
                              <right-bracket>)))))

(defn read-string [input]
  (c/matches-seq
   (h/make-state input)
   <bracket-expr>
   :success-fn (fn [product position]
                 product)
   :failure-fn (fn [error]
                 (except/throwf "SMILES parsing error: %s"
                                (h/format-parse-error error)))))
