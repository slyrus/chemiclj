
(ns chemiclj.smiles2
  (:use [chemiclj.core])
  (:require [shortcut.graph :as graph]
            [edu.arizona.fnparse [hound :as h] [core :as c]]
            [clojure.contrib [except :as except]])
  (:refer-clojure :exclude #{read-string}))

(defrecord BracketAtom
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

(h/defrule <chlorine> (h/lex (h/cat (h/lit \C) (h/lit \l))))
(h/defrule <bromine> (h/lex (h/cat (h/lit \B) (h/lit \r))))
(h/defrule <boron> (h/lex (h/cat (h/lit \B))))
(h/defrule <carbon> (h/lex (h/cat (h/lit \C))))
(h/defrule <nitrogen> (h/lex (h/cat (h/lit \N))))
(h/defrule <oxygen> (h/lex (h/cat (h/lit \O))))
(h/defrule <sulfur> (h/lex (h/cat (h/lit \S))))
(h/defrule <phosphorus> (h/lex (h/cat (h/lit \P))))
(h/defrule <fluorine> (h/lex (h/cat (h/lit \F))))
(h/defrule <iodine> (h/lex (h/cat (h/lit \I))))

(h/defrule <organic-subset-atom>
  (h/+ <chlorine> <bromine> <boron> <carbon> <nitrogen> <oxygen> <sulfur>
       <phosphorus> <fluorine> <iodine>))

(h/defrule <element-symbol>
  (h/label "an element symbol"
           (h/+ <aliphatic-element-symbol> <aromatic-element-symbol>)))

(h/defrule <atom>
  (h/+ <organic-subset-atom> <bracket-expr>))

(h/defmaker radix-natural-number [core]
  (h/hooked-rep #(+ (* core %1) %2) 0 (h/radix-digit core)))

(h/defrule <decimal-natural-number>
  (radix-natural-number 10))

(h/defrule <isotope>
  (h/label "an isotope" <decimal-natural-number>))

(h/defrule <hydrogen-count>
  (h/hook (fn [[_ hydrogen-count]]
            {:hydrogens (or hydrogen-count 1)})
          (h/cat (h/lit \H) (h/opt <decimal-natural-number>))))

(h/defrule <chirality>
  (h/hook (fn [l] {:chirality (str* (concat l))})
                    (h/+ (h/cat (h/lit \@)) (h/lit \@)
                         (h/lit \@))))

(h/defrule <charge>
  (h/+
   (h/hook (fn [[_ c]] {:charge c})
           (h/cat (h/lit \+)
                  (h/+ (h/hook (fn [_] 2) (h/lit \+))
                       (h/hook (fn [x] (or x 1))
                               (h/opt <decimal-natural-number>)))))
   (h/hook (fn [[_ c]] {:charge c})
           (h/cat (h/lit \-)
                  (h/+ (h/hook (fn [_] -2) (h/lit \-))
                       (h/hook (fn [x] (- (or x 1)))
                               (h/opt <decimal-natural-number>)))))))
(h/defrule <bracket-mods>
  (h/hook (fn [x] (when (seq x) (apply into {} (vector x))))
          (h/rep*
           (h/+
            <hydrogen-count>
            <chirality>
            <charge>))))

(h/defrule <bracket-expr>
  (h/hook (fn [[isotope symbol {:keys #{hydrogens chirality charge}} class]]
            (chemiclj.smiles2.BracketAtom.
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
   <atom>
   :success-fn (fn [product position]
                 product)
   :failure-fn (fn [error]
                 (except/throwf "SMILES parsing error: %s"
                                (h/format-parse-error error)))))
