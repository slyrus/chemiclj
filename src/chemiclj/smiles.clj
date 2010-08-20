;;; file: chemiclj/smiles.clj
;;;
;;; Copyright (c) 2010 Cyrus Harmon (ch-lisp@bobobeach.com) All rights
;;; reserved.
;;;
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
;;;
;;;   * Redistributions of source code must retain the above copyright
;;;     notice, this list of conditions and the following disclaimer.
;;;
;;;   * Redistributions in binary form must reproduce the above
;;;     copyright notice, this list of conditions and the following
;;;     disclaimer in the documentation and/or other materials
;;;     provided with the distribution.
;;;
;;; THIS SOFTWARE IS PROVIDED BY THE AUTHOR 'AS IS' AND ANY EXPRESSED
;;; OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
;;; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;;; ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
;;; DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
;;; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
;;; GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;;; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
;;; WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

(ns chemiclj.smiles
  (:use [chemiclj.core])
  (:require [chemiclj.element :as element]
            [shortcut.graph :as graph]
            [edu.arizona.fnparse [hound :as h] [core :as c]]
            [clojure.string :as str]
            [clojure.contrib [except :as except]]
            [clojure.contrib.lazy-seqs :as lazy-seqs]
            [clojure.contrib.def :as def]))

(defmacro dprint [form]
  `(let [res# ~form]
     (print res#)
     res#))

(defrecord SMILESContext
  [molecule
   last-atom
   atom-counts
   pending-rings
   open-rings
   order
   aromatic
   aromatic-atoms
   direction
   configurations])

(defn get-atom-count [context element]
  (or (get (:atom-counts context) element) 0))

(defn inc-atom-count [context element]
  (update-in context [:atom-counts element]
             #(inc (or % 0))))

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
  (h/for [symbol (h/+ (h/hook (comp str* concat) (h/cat (h/lex (h/lit \s)) (h/lit \e)))
                      (h/hook (comp str* concat) (h/cat (h/lit \a) (h/lit \s)))
                      (h/hook str (h/lit \c))
                      (h/hook str (h/lit \n))
                      (h/hook str (h/lit \o))
                      (h/hook str (h/lit \p))
                      (h/hook str (h/lit \s)))
          _ (h/alter-context
             (fn [context] (assoc context :aromatic true :hybridization :sp2)))]
         symbol))

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
          atom (h/+
                (h/for [symbol <aliphatic-organic>]
                       (let [symbol (str/capitalize symbol)]
                         (make-atom
                          symbol
                          (str symbol (inc (get-atom-count context symbol))))))
                (h/for [symbol <aromatic-organic>]
                       (let [symbol (str/capitalize symbol)]
                         (make-atom
                          symbol
                          (str symbol (inc (get-atom-count context symbol)))
                          :aromatic true
                          :hybridization :sp2))))]
         atom))

(h/defrule <bond>
  (h/+
   (h/for [_ (h/alter-context (fn [context] (assoc context :order 1)))
           bond-symbol (h/lit \-)]
          bond-symbol)
   (h/for [_ (h/alter-context (fn [context] (assoc context :order 2)))
           bond-symbol (h/lit \=)]
          bond-symbol)
   (h/for [_ (h/alter-context (fn [context] (assoc context :order 3)))
           bond-symbol (h/lit \#)]
          bond-symbol)
   (h/for [_ (h/alter-context (fn [context] (assoc context :order 4)))
           bond-symbol (h/lit \$)]
          bond-symbol)
   (h/for [_ (h/alter-context (fn [context] (assoc context :aromatic true)))
           bond-symbol (h/lit \:)]
          bond-symbol)
   (h/for [_ (h/alter-context (fn [context] (assoc context :direction :up)))
           bond-symbol (h/lit \/)]
          bond-symbol)
   (h/for [_ (h/alter-context (fn [context] (assoc context :direction :down)))
           bond-symbol (h/lit \\)]
          bond-symbol)))

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

(h/defrule <orientation>
  (h/hook (fn [config]
            {:orientation config})
          (h/+ (h/chook :clockwise
                        (h/cat (h/lex (h/lit \@)) (h/lit \@)))
               (h/chook :counter-clockwise
                        (h/lit \@)))))

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
            <orientation>
            <charge>))))

(h/defrule <bracket-expr>
  (h/for [[atom configuration]
          (h/hook (fn [[isotope symbol
                        {:keys #{hydrogens orientation charge}}
                        class context]]
                    (let [symbol (str/capitalize symbol)
                          atom (make-atom
                                symbol
                                (str symbol (inc (get-atom-count context symbol)))
                                :isotope isotope
                                :charge (or charge 0)
                                :aromatic (:aromatic context)
                                :hybridization (:hybridization context)
                                :explicit-hydrogen-count (or hydrogens 0))]
                      [atom
                       (if orientation
                         (make-tetrahedral-atom-configuration
                          atom orientation (:last-atom context) nil nil nil))]))
                  (h/label "a bracket expression"
                           (h/circumfix <left-bracket>
                                        (h/cat
                                         (h/opt <isotope>)
                                         <bracket-element-symbol>
                                         <bracket-mods>
                                         (h/opt
                                          (h/prefix (h/lit \:)
                                                    <decimal-natural-number>))
                                         h/<fetch-context>)
                                        <right-bracket>)))
          _ (h/alter-context
             (fn [context]
               (let [context (assoc context
                               :aromatic nil
                               :aromatic-atoms (conj (:aromatic-atoms context) atom))]
                 (if configuration
                   (update-in context [:configurations atom]
                              #(if %1 (conj %1 configuration) (list configuration)))
                   context))))]
         atom))

(defn- add-atom-and-bond [context atom last-atom order direction]
  (let [mol (:molecule context)]
    (let [mol
          (cond (= direction :up)
                (let [config (make-relative-vertical-configuration atom last-atom)]
                  (add-configuration
                   (add-configuration
                    mol
                    (first {atom config}))
                   (first {last-atom config})))
                (= direction :down)
                (let [config (make-relative-vertical-configuration last-atom atom)]
                  (add-configuration
                   (add-configuration
                    mol
                    (first {atom config}))
                   (first {last-atom config})))
                true mol)]
      (let [{:keys [molecule aromatic]} context]
        (cond
         (= order 1) (add-single-bond
                      (add-atom molecule atom)
                      atom last-atom)
         (= order 2) (add-double-bond
                      (add-atom molecule atom)
                      atom last-atom)
         (= order 3) (add-triple-bond
                      (add-atom molecule atom)
                      atom last-atom)
         (= order 4) (add-quadruple-bond
                      (add-atom molecule atom)
                      atom last-atom)
         true (add-bond
               (add-atom mol atom)
               atom last-atom))))))

(defn fixup-configuration [context atom last-atom]
  (if last-atom 
    (let [configuration (first
                         (filter #(= (class %1)
                                     chemiclj.core.TetrahedralAtomConfiguration)
                                 (get (:configurations context) last-atom)))]
      (if configuration
        (update-in context [:configurations last-atom]
                   (fn [x] (conj (remove #{configuration} x)
                                 (add-tetrahedral-configuration-atom configuration atom))))
        context))
    context))

(defn context-add-1-hydrogen [context atom]
  (let [mol (:molecule context)
        hatom (make-atom "H" (str "H" (inc (get-atom-count context "H"))))]
    (let [context (inc-atom-count
                   (assoc context :molecule
                          (add-bond (add-atom mol hatom) atom hatom))
                   "H")]
      (fixup-configuration context hatom atom))))

(defn context-add-n-hydrogens [context atom n]
  (loop [context context num n]
    (if (pos? num)
      (recur (context-add-1-hydrogen context atom) (dec num))
      context)))

(defn update-explicit-hydrogens [context atom]
  (let [explicit-hydrogen-count (:explicit-hydrogen-count atom)]
    (if explicit-hydrogen-count
      (context-add-n-hydrogens context atom explicit-hydrogen-count)
      context)))

(defn update-configuration [{:keys [last-atom] :as context} atom]
  (fixup-configuration context atom last-atom))

;; normally, we connect a new atom to the previous atoms
;; configuration, if any. But if we're closing a ring we also need to
;; add the current atom to the configurations of the ring opening
;; atom. Fortunately, that happens to be the new atom and the atom
;; we're connecting it to is the last-atom, so we can just reverse the
;; arguments to fixup-configuration and the right thing will happen.
;;
;; Unless of course it was the ring opening that was chiral... hmm...
(defn update-ring-atom-configuration [{:keys [last-atom] :as context} atom]
  (fixup-configuration context atom last-atom))

(defn update-last-atom [context atom]
  (assoc context :last-atom atom))

(defn post-process-atom [{:keys [last-atom] :as context} atom]
  (update-last-atom
   (update-explicit-hydrogens
    (update-configuration context atom)
    atom)
   atom))

(h/defrule <atom>
  (h/for "an atom"
         [atom (h/+ <organic-subset-atom> <bracket-expr>)
          _ (h/alter-context
             (fn [{:keys [last-atom order aromatic direction] :as context} atom]
               (assoc (inc-atom-count context (-> atom :element :id))
                 :molecule
                 (if last-atom
                   (add-atom-and-bond context atom last-atom order direction)
                   (add-atom (:molecule context) atom))
                 :order nil
                 :direction nil
                 :aromatic nil))
             atom)
          _ (h/alter-context post-process-atom atom)]
         _))

(h/defrule <ws?>
  "Consumes optional, ignored whitespace."
  (h/rep* (h/set-term "whitespace" " \t\n\r")))

(declare <chain>)

(h/defrule <branch>
  (h/for "a branch"
         [{:keys [last-atom order]} h/<fetch-context>
          branch (h/cat
                  (h/circumfix (h/lit \()
                               (h/cat
                                (h/opt (h/+ <bond> <dot>))
                                <chain>)
                               (h/lit \))))
          _ (h/alter-context
             (fn [context]
               (assoc context
                 :last-atom last-atom
                 :order order)))]
         branch))

(h/defrule <bond-or-dot>
  (h/label "a bond or dot"
           (h/+ <bond>
                <dot>)))

(defn- add-ring-bond [mol atom last-atom order]
  (cond
   (= order 1) (add-single-bond mol atom last-atom)
   (= order 2) (add-double-bond mol atom last-atom)
   (= order 3) (add-triple-bond mol atom last-atom)
   (= order 4) (add-quadruple-bond mol atom last-atom)
   true (add-bond mol atom last-atom)))

(defn- bond-symbol-order [symbol]
  (get {\- 1 \= 2 \# 3 \$ 4} symbol))

(defn- process-ring [context ring bond-symbol]
  (let [pending (get (:pending-rings context) ring)
        mol (:molecule context)
        last-atom (:last-atom context)]
    (if pending
      (let [{:keys #{atom order}} pending
            specified-order (bond-symbol-order bond-symbol)]
        (if (and (and order specified-order)
                 (not (= order specified-order)))
          (except/throwf "SMILES parsing error: %d and %d mismatch for ring bond order"
                         order specified-order))
        (let [mol (add-ring-bond mol atom last-atom (or order specified-order))]
          [(let [configuration (first
                                (filter #(= (class %1)
                                            chemiclj.core.TetrahedralAtomConfiguration)
                                        (get (:configurations context) atom)))]
             (if configuration
               (do
                 (update-in
                  context [:configurations atom]
                  (fn [x] (conj (remove #{configuration} x)
                                (replace-tetrahedral-configuration-atom
                                 configuration
                                 ring last-atom)))))
               context))
           mol (dissoc (:pending-rings context) ring) atom]))
      (do
        [(let [configuration (first
                              (filter #(= (class %1)
                                          chemiclj.core.TetrahedralAtomConfiguration)
                                      (get (:configurations context) last-atom)))]
           (if configuration
             (update-in
              context [:configurations last-atom]
              (fn [x] (conj (remove #{configuration} x)
                            (add-tetrahedral-configuration-atom configuration ring))))
             context))
         mol
         (conj (:pending-rings context)
               {ring {:atom last-atom
                      :order (bond-symbol-order bond-symbol)}})
         nil]))))

;;; FIXME Now this is broken. We should not call
;;; update-ring-atom-configuration for ring openings!!! Well, no,
;;; that's not really the problem. The problem is what to do with
;;; chiral ring opening bonds??
;;;
;;; FIXME: this context/hook stuff is a big giant mess. please clean
;;; this up.
(h/defrule <ringbond>
  (h/label "a ring bond"
           (h/for [context h/<fetch-context>
                   [context2 mol pending atom]
                   (h/hook
                    (fn [[ring-num bond-symbol]]
                      (process-ring context ring-num bond-symbol))
                    (h/+
                     (h/hook (fn [[bond _ digit1 digit2]]
                               (when (and digit1 digit2)
                                 [(+ (* 10 digit1) digit2) bond]))
                             (h/cat
                              (h/lex (h/opt <bond>)) (h/lit \%)
                              <decimal-digit> <decimal-digit>))
                     (h/hook (fn [[bond digit :as x]]
                               [digit bond])
                             (h/cat
                              (h/lex (h/opt <bond>))
                              <decimal-digit>))))
                   context (h/alter-context
                            (fn [context]
                              (assoc (update-ring-atom-configuration context2 atom)
                                :molecule mol
                                :pending-rings pending)))]
                  context)))

(h/defrule <atom-expr>
  (h/label "an atom expression"
           (h/for
            [atom <atom>
             ;; NOTE: this next rule is a hack to get around the fact
             ;; we might have a branch followed by a ring-closing at
             ;; the end of a molecule. Not sure this is allowed, but
             ;; it exists in the wild.
             _ (h/rep* <branch>)
             _ (h/rep*
                (h/for
                 [context h/<fetch-context>
                  ringbond <ringbond>]
                 ringbond))
             _ (h/alter-context
                (fn [context] (dissoc context :order)))
             _ (h/rep* <branch>)
             _ (h/opt <bond-or-dot>)
             _ (h/opt <chain>)]
            nil)))

(h/defrule <chain>
  "The main rule for a (possibly branched) chain of atoms."
  (h/label "a chain"
           (h/rep <atom-expr>)))

(defn atom-bond-orders [mol atom]
  (let [bvec (bonds mol atom)]
    (map :order bvec)))

(defn atom-aromatic-bonds [mol atom]
  (let [bvec (bonds mol atom)]
    (when (= (:hybridization atom) :sp2)
      (filter #(= (:hybridization (first (neighbors % atom))) :sp2) bvec))))

(defn atom-non-aromatic-bonds [mol atom]
  (let [bvec (bonds mol atom)]
    (when (not (= (:hybridization atom) :sp2))
      (filter #(not (= (:hybridization (first (neighbors % atom))) :sp2)) bvec))))

(defn- available-valence [mol atom]
  (- (-> atom :element element/get-normal-valences first)
     (count (neighbors mol atom))))

(defn- calculate-aromatic-bond-max-valence [mol bond]
  (let [[atom1 atom2] (atoms bond)]
    (if (and (> 2 (apply max (into (map #(or (:order %) 1) (bonds mol atom1))
                                   (map #(or (:order %) 1) (bonds mol atom2)))))
             (pos? (available-valence mol atom1))
             (pos? (available-valence mol atom2)))
      2
      1)))

(defn- fixup-sp2-atom-bonds [mol]
  (let [sp2-atoms (sort-by #(-> % :element chemiclj.element/get-normal-valences first)
                           (filter #(and
                                     (= (:hybridization %) :sp2)
                                     (seq (filter nil? (map :order (atom-aromatic-bonds mol %)))))
                                   (atoms mol)))]
    (if (seq sp2-atoms)
      (fixup-sp2-atom-bonds
       (loop [mol mol atom (first sp2-atoms)]
         (let [bondseq (filter #(nil? (:order %)) (atom-aromatic-bonds mol atom))]
           (if (seq bondseq)
             (let [mol
                   (let [atom-neighbors (neighbors mol atom)]
                     (let [valence (-> atom :element element/get-normal-valences first)]
                       (cond (> (count atom-neighbors) 3)
                             (except/throwf "Too many neighbors of sp2 atom: %s" atom)
                             
                             (< (count atom-neighbors) 2)
                             (except/throwf "Too few neighbors of sp2 atom: %s" atom)

                             (and (< (count atom-neighbors) valence)
                                  (pos? (count bondseq)))
                             (let [order (calculate-aromatic-bond-max-valence mol (first bondseq))]
                               (add-bond (remove-bond mol (first bondseq))
                                         (assoc (first bondseq)
                                           :order order)))

                             true mol)))]
               (recur mol (first (neighbors (first bondseq) atom))))
             mol))))
      mol)))

(defn fixup-non-aromatic-bonds [mol]
  (reduce #(add-bond (remove-bond %1 %2)
                     (assoc %2 :order 1))
          mol
          (filter #(nil? (:order %)) (bonds mol))))

(defn- post-process-molecule [mol]
  (-> mol
      fixup-sp2-atom-bonds
      fixup-non-aromatic-bonds))

(defn context-add-hydrogens-for-atom [context atom]
  (let [mol (:molecule context)
        valence (-> atom :element element/get-normal-valences first)
        bonds (reduce + (map :order (bonds mol atom)))
        explicit-hydrogen-count (:explicit-hydrogen-count atom)]
    (if (or (not valence) explicit-hydrogen-count)
      context
      (context-add-n-hydrogens context atom (- valence bonds)))))

(defn add-hydrogens [context]
  (let [atoms (atoms (:molecule context))]
    (reduce (fn [context atom]
              (context-add-hydrogens-for-atom context atom))
            context atoms)))

(defn read-smiles-string [input]
  (h/match
   (h/make-state input
                 :context (SMILESContext. (make-molecule) nil nil {} nil nil nil nil nil {}))
   (h/for
    [chain <chain>
     _ <ws?>
     _ h/<end-of-input>
     context h/<fetch-context>]
    context)
   :success-fn (fn [product position]
                 (reduce (fn [mol configuration]
                           (add-configuration mol configuration))
                         (fixup-non-aromatic-bonds
                          (:molecule
                           (add-hydrogens
                            (assoc product :molecule (post-process-molecule (:molecule product))))))
                         (:configurations product)))
   :failure-fn (fn [error]
                 (except/throwf "SMILES parsing error: %s"
                                (h/format-parse-error error)))))

;;; to compute the canonical SMILES we're going to need to do a few
;;; things:
;;; 1. compute the invariants for each atom in the molecule
;;; 2. assign a rank order to each atom
;;; 3. convert the rank into the nth prime
;;; 4. compute the product of the neighboring primes
;;; 5. rank the product of the primes using the previous ranks to
;;;    break ties

;;; I _think_ the SMILES spec requires us to use 1 for negative
;;; charges and 0 otherwise.
(defn charge-sgn [num] (if (neg? num) 1 0))

(defn rank [coll]
  (let [sorted (zipmap (sort (set coll)) (iterate inc 0))]
    (map #(get sorted %) coll)))

(defn rank-by [keyfn coll]
  (let [sorted (zipmap
                (sort (set (map keyfn coll)))
                (iterate inc 0))]
    (map #(get sorted (keyfn %)) coll)))

(def/defn-memo nth-prime [n]
  (nth clojure.contrib.lazy-seqs/primes n))

(defn fixpoint [coll]
  (letfn [(fixpoint* [coll last]
                     (when-let [s (seq coll)]
                       (if (= (first s) last)
                         last
                         (fixpoint* (rest s) (first s)))))]
    (fixpoint* (rest coll) (first coll))))

(defn smiles-atomic-invariant [full-molecule h-removed-molecule atom]
  (let [connections (count (graph/neighbors h-removed-molecule atom))
        non-h-bonds (reduce + (map :order (bonds h-removed-molecule atom)))
        atomic-number (:atomic-number (:element atom))
        sign-of-charge (charge-sgn (:charge atom))
        abs-charge (if (neg? (:charge atom)) (- (:charge atom)) (:charge atom))
        hydrogens (count (filter
                          #{(element/get-element "H")}
                          (map :element (graph/neighbors full-molecule atom))))]
    (+ (* 10000000 connections)
       (* 100000 non-h-bonds)
       (* 1000 atomic-number)
       (* 100 sign-of-charge)
       (* 10 abs-charge)
       hydrogens)))

(defn smiles-atomic-invariants [full-mol]
  (let [mol (remove-atoms-of-element full-mol "H")]
    (reduce (fn [m atom]
              (assoc m atom (smiles-atomic-invariant full-mol mol atom)))
            {}
            (atoms mol))))

(defn smiles-atomic-invariant-ranks [mol]
  (let [invariants (smiles-atomic-invariants mol)]
    (zipmap
     (map first invariants)
     (map inc (rank-by second invariants)))))

;; note that ranks are 0-indexed, but the ranks in the weininger paper
;; start at 1. shouldn't make a difference in the answer, but if one
;; inspects the ranks, beware of this. Note that we get the same
;; multiple of primes either way.
(defn ranks-and-prime-products [mol imap]
  (reduce (fn [m [atom rank]]
            (assoc m atom
                   [rank (reduce * (map #(or (when-let [p (get imap %)]
                                               (nth-prime (dec p)))
                                             1)
                                        (graph/neighbors mol atom)))]))
          {}
          imap))

(defn rank-coll-by-second [coll]
  (zipmap
   (map first coll)
   (map inc (rank-by second coll))))

;;; FIXME! We need to fix this such that we don't sort the atoms by
;;; name, which means we need to preserve the order of atoms between
;;; iterations of the canonical labeling algorithm.
(defn break-ties [ranked-atoms]
  (let [freqs (frequencies (map second ranked-atoms))]
    (let [lowest (ffirst (filter #(> (val %) 1) freqs))]
      (if lowest
        (first (reduce (fn [[acc lowest] [atom rank]]
                         (if (= rank lowest)
                           [(conj acc {atom (dec (* 2 rank))}) nil]
                           [(conj acc {atom (* 2 rank)}) lowest]))
                       [{} lowest]
                       ranked-atoms))
        ranked-atoms))))

(defn smiles-canonical-labels [mol]
  (let [ranks (smiles-atomic-invariant-ranks mol)]
    (fixpoint
     (iterate (fn [ranks]
                (let [sums (ranks-and-prime-products mol ranks)]
                  (break-ties (rank-coll-by-second sums))))
              ranks))))

;;; TODO

(def *organic-subset-elements*
     (set (map element/get-element ["B" "C" "N" "O" "P" "S" "F" "Cl" "Br" "I"])))

(defn organic-subset? [atom]
  (*organic-subset-elements* atom))

;; the problem here is that we need to use different criteria for
;; ordering the neighbors if we are in a ring, or not. Or perhaps if
;; we're about to close a ring. I'm not sure which... The smiles paper
;; is unclear on this.
(defn- smiles-neighbors [mol labels atom]
  (let [bonds (bonds mol atom)]
    (map first (map #(graph/neighbors % atom)
                    (sort-by #(vector ((comp - :order) %)
                                      (labels (first (graph/neighbors % atom))))
                             bonds)))))

(defn first-pass [mol labels atom visited rings]
  (let [visited (conj visited atom)]
    (loop [neighbors (filter (complement visited)
                             (smiles-neighbors mol labels atom))
           mol mol visited visited rings rings]
      (if (seq neighbors)
        (let [neighbor (first neighbors)]
          (let [[mol visited rings]
                (first-pass
                 (remove-bond mol atom neighbor)
                 labels
                 neighbor
                 visited
                 (if (visited neighbor)
                   (conj rings [atom neighbor])
                   rings))]
            (recur (rest neighbors) mol visited rings)))
        [mol visited rings]))))

(defn first-empty [set]
  (let [x (first (filter (fn [[a b]] (not (= a b)))
                         (map vector
                              (iterate inc 1)
                              (sort (keys set)))))]
    (or (first x) (inc (count set)))))

(defn write-bond [bond]
  (when bond
    (cond (= (:order bond) 2)
          (print "="))))

(defn write-ring-number [number]
  (if (< number 10)
    (print number)
    (print "%" number)))

(defn write-bracket-atom [atom]
  (let [element (:element atom)
        hydrogens (:explicit-hydrogen-count atom)
        charge (:charge atom)
        isotope (:isotope atom)]
    (print "[")
    (when isotope
      (print isotope))
    (print (:id element))
    (cond (= hydrogens 1)
          (print "H")
          (> hydrogens 1)
          (do
            (print "H")
            (print hydrogens)))
    (cond (= 1 charge)
          (print "+")
          (pos? charge)
          (dotimes [i charge]
            (print "+"))
          (= -1 charge)
          (print "-")
          (neg? charge)
          (dotimes [i (- charge)]
            (print "-")))
    (print "]")))

(defn write-atom [atom]
  (let [element (:element atom)
        hydrogens (:explicit-hydrogen-count atom)
        charge (:charge atom)
        isotope (:isotope atom)]
    (if (and (organic-subset? element)
             (nil? hydrogens)
             (or (nil? charge) (zero? charge))
             (nil? isotope))
      (print (:id element))
      (write-bracket-atom atom))))

(def *reuse-ring-number* false)

(defn second-pass [mol labels atom bond visited rings open-rings ring-count]
  (if (visited atom)
    [mol visited rings open-rings]
    (let [visited (conj visited atom)]
      (write-bond bond)
      (write-atom atom)
      (let [[open-rings ring-bonds ring-count]
            (reduce (fn [[open-rings ring-bonds ring-count] ring]
                      (let [ring-num (if *reuse-ring-number*
                                       (first-empty open-rings)
                                       (inc ring-count))]
                        (write-ring-number ring-num)
                        [(assoc open-rings ring-num ring)
                         (conj ring-bonds (bond? mol (first ring) (second ring)))
                         (inc ring-count)]))
                    [open-rings #{} ring-count]
                    (filter #(= (first %) atom) rings))]
        ;; now for the ring closings
        (let [open-rings
              (reduce (fn [open-rings [ring-num ring]]
                        (let [ring-bond (bond? mol (first ring) (second ring))]
                          (write-bond ring-bond))
                        (write-ring-number ring-num)
                        (dissoc open-rings ring-num))
                      open-rings
                      (sort-by key (filter #(= (second (val %)) atom) open-rings)))]
          (loop [neighbors (filter (complement visited)
                                   (smiles-neighbors (reduce (fn [mol bond]
                                                               (remove-bond mol bond))
                                                             mol
                                                             ring-bonds)
                                                     labels
                                                     atom))
                 mol mol visited visited rings rings open-rings open-rings ring-count ring-count]
            (if (seq neighbors)
              (let [branch (seq (rest neighbors))]
                (when branch (print "("))
                (let [neighbor (first neighbors)]
                  (let [[mol visited rings open-rings ring-count]
                        (second-pass
                         (remove-bond mol atom neighbor)
                         labels
                         neighbor
                         (bond? mol atom neighbor)
                         visited
                         (if (visited neighbor)
                           (conj rings [neighbor atom])
                           rings)
                         open-rings
                         ring-count)]
                    (when branch (print ")"))
                    (recur (rest neighbors) mol visited rings open-rings ring-count))))
              [mol visited rings open-rings ring-count])))))))

(defn compute-explicit-hydrogens [mol]
  ;; TODO add explicit hydrogen count here as necessary!
  mol)

(defn remove-simple-hydrogens [mol]
  (reduce (fn [mol atom]
            (if (:isotope atom)
              mol
              (remove-atom mol atom)))
          mol
          (get-atoms-of-element mol "H")))

(defn write-smiles-string [molecule]
  ;; TODO special cases for H and H2
  (let [mol (remove-simple-hydrogens
             (compute-explicit-hydrogens molecule))
        labels (smiles-canonical-labels mol)
        start (ffirst (sort-by second labels))]
    (let [[new-mol _ rings] (first-pass mol labels start #{} #{})]
      (with-out-str
        (second-pass mol labels start nil #{} rings {} 0)))))
