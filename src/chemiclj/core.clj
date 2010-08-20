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
  (:require [chemiclj.element :as element]
            [shortcut.graph :as g]
            [clojure.string :as string]
            [clojure.contrib.def :as d]
            [clojure.contrib [except :as except]]))

(d/defalias neighbors g/neighbors)

(defprotocol PAtomContainer
  (atoms [obj]))

(defprotocol PMass
  (mass [mol])
  (exact-mass [mol]))

(defprotocol PMolecule
  (add-atom [mol atom])
  (remove-atom [mol atom])
  (bonds [mol] [mol atom])
  (bond? [mol atom1 atom2])
  (add-bond [mol bond] [mol atom1 atom2])
  (remove-bond [mol bond] [mol atom1 atom2])
  (configurations [mol])
  (add-configuration [mol configuration]))


;;;
;;; Atom record and related functions
(defrecord Atom [_name element isotope chirality charge
                 hybridization aromatic explicit-hydrogen-count]
  PMass
  (mass [atm] (-> atm :element :mass))
  (exact-mass [atm]
              (if (:isotope atm)
                (:exact-mass (:isotope atm))
                (:exact-mass (first (element/element-abundnant-isotopes (:element atm))))))

  clojure.lang.Named
  (getName [atm] _name))

(defn make-atom [element name & {:keys [isotope chirality charge
                                        hybridization aromatic
                                        explicit-hydrogen-count]
                                 :or {charge 0}}]
  (Atom. name (element/get-element element)
         isotope chirality charge hybridization
         aromatic explicit-hydrogen-count))

(defn get-atom [mol atom]
  (cond (= (type atom) java.lang.String)
        (when-let [res (filter #(= (name %) atom) (atoms mol))]
          (first res))
        true
        (get (atoms mol) atom)))


;;;
;;; Bond record and related functions
(defrecord Bond [_nodes type order direction]
  PAtomContainer
  (atoms [bond] _nodes)

  g/NodeSet
  (g/nodes [bond] _nodes)
  (g/node? [bond node] (some #{node} _nodes))
  (g/neighbors [bond node] (remove #{node} _nodes))

  g/Edge
  (g/left [bond] (first _nodes))
  (g/right [bond] (second _nodes)))

(defn make-bond [atom1 atom2 & {:keys [type order direction]}]
  (Bond. [atom1 atom2] type order direction))

;;;
;;; Molecule record and related functions
(defrecord Molecule [_graph _name _configurations]
  PAtomContainer
  (atoms [mol] (g/nodes _graph))

  PMolecule
  (add-atom [mol atom]
            (assoc mol :_graph (g/add-node _graph atom)))
  (remove-atom [mol atom]
               (assoc mol :_graph (g/remove-node _graph atom)))
  (bonds [mol] (g/edges _graph))
  (bonds [mol atom] (g/edges _graph (get-atom mol atom)))
  (bond? [mol atom1 atom2] (g/edge? _graph atom1 atom2))
  (add-bond [mol bond]
            (assoc mol :_graph (g/add-edge _graph bond)))
  (add-bond [mol atom1 atom2]
            (assoc mol :_graph
                   (g/add-edge _graph
                               (apply make-bond
                                      (get-atom mol atom1)
                                      (get-atom mol atom2)
                                      nil))))
  (remove-bond [mol bond]
               (assoc mol :_graph (g/remove-edge _graph bond)))
  (remove-bond [mol atom1 atom2]
               (assoc mol :_graph (g/remove-edge _graph atom1 atom2)))
  (configurations [mol] _configurations)
  (add-configuration [mol configuration]
                     (update-in mol [:_configurations (key configuration)]
                                #(into % (val configuration))))

  PMass
  (mass [mol] (reduce + (map mass (atoms mol))))
  (exact-mass [mol] (reduce + (map exact-mass (atoms mol))))

  g/NodeSet
  (g/nodes [mol] (g/nodes _graph))
  (g/node? [mol node] (g/node? _graph node))
  (g/neighbors [mol node] (g/neighbors _graph node))

  clojure.lang.Named
  (getName [mol] _name))

(defn count-elements [mol]
  (reduce
   (fn [counts atom]
     (assoc counts (:element atom) (inc (or (get counts (:element atom)) 0))))
   {}
   (shortcut.graph/breadth-first-traversal (:_graph mol) (first (atoms mol)))))

(defn molecular-formula [mol]
  (apply str
         (map #(str (:id (first %)) (second %))
              (sort-by #(:id (first %)) (chemiclj.core/count-elements mol)))))

(defn make-molecule
  ([]
     (Molecule. (g/make-graph) nil {}))
  ([atoms atom-pairs-vec]
     (make-molecule atoms atom-pairs-vec nil))
  ([atoms atom-pairs-vec name]
     (reduce
      (fn [mol [atom1 atom2]]
        (assoc mol :_graph (g/add-edge (:_graph mol) atom1 atom2 (make-bond atom1 atom2))))
      (Molecule. (g/make-graph atoms) nil {})
      atom-pairs-vec)))

(defn name-molecule [mol name]
  (conj mol {:_name name}))

(defn names [seq]
  (map name seq))

;;; should these be protocol methods???
(defn add-single-bond [mol atom1 atom2]
  (add-bond mol (make-bond atom1 atom2 :type :single :order 1)))

(defn add-double-bond [mol atom1 atom2]
  (add-bond mol (make-bond atom1 atom2 :type :double :order 2)))

(defn add-triple-bond [mol atom1 atom2]
  (add-bond mol (make-bond atom1 atom2 :type :triple :order 3)))

(defn add-quadruple-bond [mol atom1 atom2]
  (add-bond mol (make-bond atom1 atom2 :type :quadruple :order 4)))

(defn add-aromatic-bond [mol atom1 atom2]
  (add-bond mol (make-bond atom1 atom2 :type :aromatic :order 1.5)))

(defn add-atom* [mol element name & [attached-to]]
  (let [atm (make-atom element name)]
    (if attached-to
      (add-bond (add-atom mol atm) (make-bond atm attached-to))
      (add-atom mol atm))))

(defrecord TetrahedralAtomConfiguration [center direction w x y z])

(defn make-tetrahedral-atom-configuration [center direction w x y z]
  (TetrahedralAtomConfiguration. center direction w x y z))

(defn get-tetrahedral-configuration-vector [configuration]
  [(:w configuration) (:x configuration) (:y configuration) (:z configuration)])

(defn set-tetrahedral-configuration-vector [configuration v]
  (assoc configuration
    :w (nth v 0)
    :x (nth v 1)
    :y (nth v 2)
    :z (nth v 3)))

(defn add-tetrahedral-configuration-atom [configuration atom]
  (cond (nil? (:w configuration)) (assoc configuration :w atom)
        (nil? (:x configuration)) (assoc configuration :x atom)
        (nil? (:y configuration)) (assoc configuration :y atom)
        (nil? (:z configuration)) (assoc configuration :z atom)
        true (except/throwf "Too many neighbors for tetrahedral atom %s"
                            atom)))

(defn replace-val [m old new]
  (reduce (fn [m [k v]]
            (if (= v old)
              (assoc m k new)
              m))
          m
          m))

(defn replace-tetrahedral-configuration-atom [configuration old new]
  (reduce (fn [m [k v]]
            (if (= v old)
              (assoc m k new)
              m))
          configuration
          configuration))

(defrecord RelativeVerticalConfiguration [top bottom])

(defn make-relative-vertical-configuration [top bottom]
  (RelativeVerticalConfiguration. top bottom))

(defrecord DoubleBondConfiguration [a b c d])

(defn get-atoms-of-element [mol elmnt]
  (let [elmnt (element/get-element elmnt)]
    (filter #(= (:element %) elmnt) (atoms mol))))

(defn remove-atoms-of-element [mol element]
  (reduce (fn [mol atom]
            (remove-atom mol atom))
          mol
          (get-atoms-of-element mol element)))

