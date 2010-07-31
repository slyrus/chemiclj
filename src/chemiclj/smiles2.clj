
(ns chemiclj.smiles2
  (:use [chemiclj.core])
  (:require [shortcut.graph :as graph]
            [edu.arizona.fnparse [hound :as h] [core :as c]]
            [clojure.string :as str]
            [clojure.contrib [except :as except]]))

(defrecord SMILESContext
  [molecule last-atom atom-counts rings])

(defn get-atom-count [context element]
  (or (get (:atom-counts context) element) 0))

(defn inc-atom-count [context element]
  (assoc context :atom-counts
         (assoc (:atom-counts context) element
                (inc (get-atom-count context element)))))

(defrecord BracketAtom
  [isotope symbol chirality hydrogen-count charge class])

(defn str* [objects]
  (apply str objects))

(h/defrule <left-bracket> (h/lit \[))
(h/defrule <right-bracket> (h/lit \]))

(h/defrule <bracket-aliphatic-element-symbol>
  (h/hook (comp str* concat)
          (h/cat h/<uppercase-ascii-letter>
                 (h/opt h/<lowercase-ascii-letter>)
                 (h/opt h/<lowercase-ascii-letter>))))

(h/defrule <bracket-aromatic-element-symbol>
  (h/+ (h/hook (comp str* concat) (h/cat (h/lit \s) (h/lit \e)))
       (h/hook (comp str* concat) (h/cat (h/lit \a) (h/lit \s)))
       (h/hook str (h/lit \c))
       (h/hook str (h/lit \n))
       (h/hook str (h/lit \o))
       (h/hook str (h/lit \p))
       (h/hook str (h/lit \s))))

(h/defrule <bracket-element-symbol>
  (h/label "a bracket element symbol"
           (h/+ <bracket-aliphatic-element-symbol> <bracket-aromatic-element-symbol>)))

(h/defrule <aliphatic-chlorine> (h/hook (comp str* concat)
                              (h/cat (h/lex (h/lit \C)) (h/lit \l))))
(h/defrule <aliphatic-bromine> (h/hook (comp str* concat)
                             (h/cat (h/lex (h/lit \B)) (h/lit \r))))
(h/defrule <aliphatic-boron> (h/hook str (h/lit \B)))
(h/defrule <aliphatic-carbon> (h/hook str (h/lit \C)))
(h/defrule <aliphatic-nitrogen> (h/hook str (h/lit \N)))
(h/defrule <aliphatic-oxygen> (h/hook str (h/lit \O)))
(h/defrule <aliphatic-sulfur> (h/hook str (h/lit \S)))
(h/defrule <aliphatic-phosphorus> (h/hook str (h/lit \P)))
(h/defrule <aliphatic-fluorine> (h/hook str (h/lit \F)))
(h/defrule <aliphatic-iodine> (h/hook str (h/lit \I)))

(h/defrule <aliphatic-organic>
  (h/+ <aliphatic-chlorine>
       <aliphatic-bromine>
       <aliphatic-boron>
       <aliphatic-carbon>
       <aliphatic-nitrogen>
       <aliphatic-oxygen>
       <aliphatic-sulfur>
       <aliphatic-phosphorus>
       <aliphatic-fluorine>
       <aliphatic-iodine>))

(h/defrule <aromatic-boron> (h/hook str (h/lit \b)))
(h/defrule <aromatic-carbon> (h/hook str (h/lit \c)))
(h/defrule <aromatic-nitrogen> (h/hook str (h/lit \n)))
(h/defrule <aromatic-oxygen> (h/hook str (h/lit \o)))
(h/defrule <aromatic-sulfur> (h/hook str (h/lit \s)))
(h/defrule <aromatic-phosphorus> (h/hook str (h/lit \p)))

(h/defrule <aromatic-organic>
  (h/+ <aromatic-boron>
       <aromatic-carbon>
       <aromatic-nitrogen>
       <aromatic-oxygen>
       <aromatic-sulfur>
       <aromatic-phosphorus>))

(h/defrule <organic-subset-atom>
  (h/for [context h/<fetch-context>
          atom (h/hook
                (fn [symbol]
                  (make-atom symbol
                             (str symbol (inc (get-atom-count context (str/upper-case symbol))))))
                (h/+ <aliphatic-organic> <aromatic-organic>))]
         atom))

(h/defrule <bond>
  (h/+
   (h/lit \-)
   (h/lit \=)
   (h/lit \#)
   (h/lit \$)
   (h/lit \:)
   (h/lit \/)
   (h/lit \\)))

(h/defrule <dot> (h/lit \.))

(h/defrule <decimal-digit> (h/radix-digit 10))

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
  (h/for [context h/<fetch-context>
          atom (h/hook (fn [[isotope symbol {:keys #{hydrogens chirality charge}} class]]
                         (make-atom symbol (str symbol (inc (get-atom-count context symbol)))
                                    isotope chirality charge nil hydrogens))
                       (h/label "a bracket expression"
                                (h/circumfix <left-bracket>
                                             (h/cat
                                              (h/opt <isotope>)
                                              <bracket-element-symbol>
                                              <bracket-mods>
                                              (h/opt
                                               (h/prefix (h/lit \:)
                                                         <decimal-natural-number>)))
                                             <right-bracket>)))]
         atom))

(h/defrule <atom>
  (h/for "an atom"
         [atom (h/+ <organic-subset-atom> <bracket-expr>)
          _ (h/alter-context
             (fn [{:keys [last-atom] :as context} atom]
               (assoc (inc-atom-count context (-> atom :element :id))
                 :molecule (if last-atom
                             (add-bond (add-atom (:molecule context) atom)
                                       atom last-atom)
                             (add-atom (:molecule context) atom))
                 :last-atom atom))
             atom)]
         atom))

(h/defrule <ws?>
  "Consumes optional, ignored whitespace."
  (h/rep* (h/set-term "whitespace" " \t\n\r")))

(h/defrule <ringbond>
  (h/label "a ringbond"
           (h/+
            (h/cat
             (h/lex (h/opt <bond>)) (h/lit \%) <decimal-digit> <decimal-digit>)
            (h/cat
             (h/lex (h/opt <bond>)) <decimal-digit>))))

(declare <branch>)
(declare <chain>)

(h/defrule <bond-or-dot>
  (h/label "a bond or a dot"
           (h/+ <bond> <dot>)))

(h/defrule <branched-atom>
  (h/label "a branched-atom"
           (h/cat <atom>
                  (h/rep* <ringbond>)
                  (h/rep* <branch>)
                  (h/opt <bond-or-dot>)
                  (h/opt <chain>))))

(declare <chain>)

(h/defrule <branch>
  (h/for "a branch"
         [{:keys [last-atom]} h/<fetch-context>
          branch (h/cat
                  (h/circumfix (h/lit \()
                               (h/cat
                                (h/opt (h/+ <bond> <dot>))
                                <chain>)
                               (h/lit \))))
          _ (h/alter-context
             (fn [context]
               (assoc context :last-atom last-atom)))]
         branch))

(h/defrule <chain>
  (h/label "a chain"
           (h/rep <branched-atom>)))

(defn read-smiles-string [input]
  (h/match
   (h/make-state input
                 :context (SMILESContext. (make-molecule) nil nil nil))
   (h/for
    [chain <chain>
     _ <ws?>
     _ h/<end-of-input>
     context h/<fetch-context>]
    context)
   :success-fn (fn [product position]
                 product)
   :failure-fn (fn [error]
                 (except/throwf "SMILES parsing error: %s"
                                (h/format-parse-error error)))))
