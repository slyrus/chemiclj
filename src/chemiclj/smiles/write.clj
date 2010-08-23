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

(ns chemiclj.smiles.write
  (:use [chemiclj core atom bond molecule])
  (:require [chemiclj.element :as element]
            [shortcut.graph :as graph]
            [clojure.contrib.lazy-seqs :as lazy-seqs]
            [clojure.contrib.def :as def]))

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
          (print "="))
    (cond (= (:order bond) 3)
          (print "#"))))

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
