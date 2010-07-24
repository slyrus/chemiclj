
(ns chemiclj.core
  (:use chemiclj.element)

  (:use shortcut.graph)

  (:require [clojure.string :as string]))

(defrecord Atom [name element charge isotope hybridization])

(defn atom-exact-mass [atom]
  (if (:isotope atom)
    (:exact-mass (:isotope atom))
    (:exact-mass (first (element-abundnant-isotopes (:element atom))))))

(defrecord Bond [nodes type order direction]
  NodeSet
  (nodes [bond] nodes)
  (node? [bond node] (some #{node} nodes))
  (neighbors [bond node] (remove #{node} nodes))

  clojure.lang.Indexed
  (nth [bond i] (nth nodes i)))

(defprotocol PMolecule
  (atoms [mol])
  (bonds [mol])
  (mass [mol])
  (exact-mass [mol]))

(defrecord Molecule [_graph _name]
  PMolecule
  (atoms [mol] (nodes _graph))
  (bonds [mol] (edges _graph))
  (mass [mol]
        (reduce + (map #(-> % :element :mass) (atoms mol))))
  (exact-mass [mol]
              (reduce + (map atom-exact-mass (atoms mol)))))

(defn make-atom
  ([element name]
     (Atom. name (get-element element) nil nil nil)))

(defn make-bond [atom1 atom2 & {:keys [type order direction],
                                :or {type :single order 1}}]
  (Bond. [atom1 atom2] type order direction))

(defn make-molecule [atoms atom-pairs-vec & name]
  (reduce
  (fn [g [atom1 atom2]]
     (add-edge g atom1 atom2 (make-bond atom1 atom2)))
  (Molecule. (make-graph atoms) name)
  atom-pairs-vec))
