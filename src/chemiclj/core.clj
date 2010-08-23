;;; file: chemiclj/core.clj
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

(ns chemiclj.core
  (:require [clojure.contrib.def :as d]
            [shortcut.graph :as g]))

(d/defalias neighbors g/neighbors)

(defprotocol PElement
  (isotopes [obj]))

(defprotocol PMass
  (mass [mol])
  (exact-mass [mol]))

(defprotocol PAtomContainer
  (atoms [obj]))

(defprotocol PMolecule
  (add-atom [mol atom])
  (remove-atom [mol atom])
  (bonds [mol] [mol atom])
  (bond? [mol atom1 atom2])
  (add-bond [mol bond] [mol atom1 atom2])
  (remove-bond [mol bond] [mol atom1 atom2])
  (configurations [mol])
  (add-configuration [mol configuration]))

(defn names [seq]
  (map name seq))

(defn get-atom [mol atom]
  (cond (= (type atom) java.lang.String)
        (when-let [res (filter #(= (name %) atom) (atoms mol))]
          (first res))
        true
        (get (atoms mol) atom)))

(use '[chemiclj.element])

(use '[chemiclj.atom])
(defn make-atom [element name & {:keys [isotope chirality charge
                                        hybridization aromatic
                                        explicit-hydrogen-count]
                                 :or {charge 0}}]
  (chemiclj.atom.Atom. name (get-element element)
                       isotope chirality charge hybridization
                       aromatic explicit-hydrogen-count))


(use '[chemiclj.bond])
(defn make-bond [atom1 atom2 & {:keys [type order direction]}]
  (chemiclj.bond.Bond. [atom1 atom2] type order direction))


(use '[chemiclj.molecule])
(defn molecular-formula [mol]
  (apply str
         (map #(str (:id (first %)) (second %))
              (sort-by #(:id (first %)) (count-elements mol)))))

(defn make-molecule
  ([]
     (chemiclj.molecule.Molecule. (g/make-graph) nil {}))
  ([atoms atom-pairs-vec]
     (make-molecule atoms atom-pairs-vec nil))
  ([atoms atom-pairs-vec name]
     (reduce
      (fn [mol [atom1 atom2]]
        (assoc mol :_graph (g/add-edge (:_graph mol) atom1 atom2 (make-bond atom1 atom2))))
      (chemiclj.molecule.Molecule. (g/make-graph atoms) nil {})
      atom-pairs-vec)))

(defn name-molecule [mol name]
  (conj mol {:_name name}))

(defn get-atoms-of-element [mol elmnt]
  (let [elmnt (get-element elmnt)]
    (filter #(= (:element %) elmnt) (atoms mol))))

(defn remove-atoms-of-element [mol element]
  (reduce (fn [mol atom]
            (remove-atom mol atom))
          mol
          (get-atoms-of-element mol element)))

