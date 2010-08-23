;;; file: test/chemiclj_scratch.clj
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


(ns chemiclj-scratch
  (:use [chemiclj.core])
  (:require [chemiclj.element :as element]
            [chemiclj.smiles :as smiles]
            [chemiclj.smiles.write :as write]
            [smiles-test :as smiles-test]
            [shortcut.graph :as graph]))

(def c1 (make-atom :c "C1"))
(def h1 (make-atom :h "H1"))
(def h2 (make-atom :h "H2"))
(def h3 (make-atom :h "H3"))
(def c2 (make-atom :c "C2"))
(def h4 (make-atom :h "H4"))
(def h5 (make-atom :h "H5"))
(def h6 (make-atom :h "H6"))

(def ethane (make-molecule #{c1 h1 h2 h3 c2 h4 h5 h6}
                           [[c1 h1]
                            [c1 h2]
                            [c1 h3]
                            [c1 c2]
                            [c2 h4]
                            [c2 h5]
                            [c2 h6]]))
(mass ethane)
(exact-mass ethane)

(add-bond (reduce add-atom (make-molecule) [c1 c2])
          (make-bond c1 c2))

(add-bond (reduce add-atom (make-molecule) [c1 c2]) c1 c2)



;;; scratch for writing SMILES strings

(let [mol (smiles-test/get-molecule "6-amino-2-ethyl-5-(aminomethyl)-1-hexanol")]
  (write/smiles-atomic-invariants mol))

(let [mol (smiles-test/get-molecule "6-amino-2-ethyl-5-(aminomethyl)-1-hexanol")]
  (write/rank-by second (write/smiles-atomic-invariants mol)))

(reduce #(assoc %1 ((comp name first) %2) (second %2))
                 {}
                 (let [mol (smiles-test/get-molecule
                            "6-amino-2-ethyl-5-(aminomethyl)-1-hexanol")]
                   (write/smiles-atomic-invariants mol)))

(reduce #(assoc %1 ((comp name first) %2) (second %2))
                 {}
                 (let [mol (smiles-test/get-molecule
                            "6-amino-2-ethyl-5-(aminomethyl)-1-hexanol")]
                   (write/smiles-atomic-invariant-ranks mol)))

;;;

(sort-by first
         (map #(vector ((comp name first) %)
                       ((comp inc second) %))
              (let [mol (smiles-test/get-molecule
                         "6-amino-2-ethyl-5-(aminomethyl)-1-hexanol")]
                (write/smiles-canonical-labels mol))))

(let [mol (smiles-test/get-molecule "6-amino-2-ethyl-5-(aminomethyl)-1-hexanol")
      labels (write/smiles-canonical-labels mol)]
  (graph/depth-first-traversal mol
                               (ffirst (sort-by second labels))))


(let [mol (smiles-test/get-molecule
           "6-amino-2-ethyl-5-(aminomethyl)-1-hexanol")
      labels (write/smiles-canonical-labels mol)]
  (ffirst (sort-by second labels)))


(sort-by first
         (map #(vector ((comp name first) %)
                       (second %))
              (let [mol (smiles-test/get-molecule
                         "cubane")]
                (write/smiles-canonical-labels mol))))

(smiles/read-smiles-string "Br/C(=C/F)I")

(map (fn [[k v]] [(name k) (map (fn [x] (names [(:top x)
                                                (:bottom x)])) v)])
     (configurations (smiles/read-smiles-string "F/C=C/F")))

(map (fn [[k v]] [(name k) (map (fn [x] (names [(:top x)
                                                (:bottom x)])) v)])
     (configurations (smiles/read-smiles-string "Br/C(=C\\F)/I")))

(= (set
    (map (fn [[k v]] [(name k) (map (fn [x] (names [(:top x)
                                                    (:bottom x)])) v)])
         (configurations (smiles/read-smiles-string "Br/C(/I)=C/F"))))

   (set
    (map (fn [[k v]] [(name k) (map (fn [x] (names [(:top x)
                                                    (:bottom x)])) v)])
         (configurations (smiles/read-smiles-string "Br/C(=C/F)/I")))))

(map (fn [[k v]] [(name k) (map (fn [x] (names [(:top x)
                                                (:bottom x)])) v)])
     (configurations (smiles-test/get-molecule "tamoxifen")))

;; find the double bonds in tamoxifen:
(filter (comp #{2} :order) (bonds (smiles-test/get-molecule "tamoxifen")))

;; get the names of all the atoms that participate in double bonds in tamoxifen:
(map (comp names atoms)
     (filter (comp #{2} :order) (bonds (smiles-test/get-molecule "tamoxifen"))))

;; or
(names
 (reduce #(into %1 (atoms %2))
         []
         (filter (comp #{2} :order) (bonds (smiles-test/get-molecule "tamoxifen")))))

;; now let's get the top atom of each configuration attached to each
;; atom that participates in a double bond:
(let [mol (smiles-test/get-molecule "tamoxifen")]
  (map (partial map (comp name :top))
       (map second
            (let [s (set (reduce #(into %1 (atoms %2))
                                 []
                                 (filter (comp #{2} :order) (bonds mol))))]
              (filter #(s (first %)) (configurations mol))))))


;;; oxytocin
(map (fn [[k v]]
       [(name k)
        (map (fn [x]
               (cond (= (class x) chemiclj.molecule.RelativeVerticalConfiguration)
                     (names [(:top x) (:bottom x)])
                     (= (class x) chemiclj.molecule.TetrahedralAtomConfiguration)
                     [(name (:center x))
                      (:direction x)
                      (names [(:w x) (:x x) (:y x) (:z x)])]))
             v)])
     (configurations (smiles-test/get-molecule "oxytocin")))

;;;
(map (fn [[k v]]
       [(name k)
        (map (fn [x]
               (cond (= (class x) chemiclj.molecule.RelativeVerticalConfiguration)
                     (names [(:top x) (:bottom x)])
                     (= (class x) chemiclj.molecule.TetrahedralAtomConfiguration)
                     [(name (:center x))
                      (:direction x)
                      (names [(:w x) (:x x) (:y x) (:z x)])]))
             v)])
     (configurations (smiles/read-smiles-string "C[C@H](Br)I")))

(map (fn [[k v]]
       [(name k)
        (map (fn [x]
               (cond (= (class x) chemiclj.molecule.RelativeVerticalConfiguration)
                     (names [(:top x) (:bottom x)])
                     (= (class x) chemiclj.molecule.TetrahedralAtomConfiguration)
                     [(name (:center x))
                      (:direction x)
                      (names [(:w x) (:x x) (:y x) (:z x)])]))
             v)])
     (configurations (smiles/read-smiles-string "N1C[C@H]1C")))

(map (fn [[k v]]
       [(name k)
        (map (fn [x]
               (cond (= (class x) chemiclj.molecule.RelativeVerticalConfiguration)
                     (names [(:top x) (:bottom x)])
                     (= (class x) chemiclj.molecule.TetrahedralAtomConfiguration)
                     [(name (:center x))
                      (:direction x)
                      (names [(:w x) (:x x) (:y x) (:z x)])]))
             v)])
     (configurations (smiles/read-smiles-string "[H][C@]1(Br)CCC1")))

;;; names of each non-H atom and its neighbors:
(let [mol (smiles-test/get-molecule "cubane")]
  (map (fn [a]
     [(name a) (names (neighbors mol a))])
       (atoms (remove-atoms-of-element mol "H"))))

(let [mol (smiles-test/get-molecule "serotonin")]
  (map (fn [a]
     [(name a) (names (neighbors mol a))])
       (atoms (remove-atoms-of-element mol "H"))))


(map (comp names atoms)
     (filter (comp #{2} :order) (bonds (smiles-test/get-molecule "tropone"))))


(map (comp names atoms) (bonds (smiles/read-smiles-string "[Na+].[O-]c1ccccc1")))
(map (comp names atoms) (bonds (smiles/read-smiles-string "c1cc([O-].[Na+])ccc1")))
