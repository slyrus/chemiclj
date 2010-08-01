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
  (:use [chemiclj.element :as element]
        [shortcut.graph
         :only (NodeSet Edge
                make-graph
                nodes add-node add-nodes
                edges add-edge add-edges
                neighbors
                connected-components connected-component partition-graph)
         :as graph])
  (:require [clojure.string :as string]))

(defprotocol PHasAtoms
  (atoms [obj]))

(defprotocol PMass
  (mass [mol])
  (exact-mass [mol]))

(defrecord Atom [_name element isotope chirality
                 charge hybridization explicit-hydrogen-count]
  PMass
  (mass [atm] (-> atm :element :mass))
  (exact-mass [atm]
              (if (:isotope atm)
                (:exact-mass (:isotope atm))
                (:exact-mass (first (element-abundnant-isotopes (:element atm))))))

  clojure.lang.Named
  (getName [atm] _name))

(defrecord Bond [_nodes type order direction]
  PHasAtoms
  (atoms [bond] _nodes)

  NodeSet
  (nodes [bond] _nodes)
  (node? [bond node] (some #{node} _nodes))
  (neighbors [bond node] (remove #{node} _nodes))

  Edge
  (left [bond] (first _nodes))
  (right [bond] (second _nodes)))

(defprotocol PMolecule
  (bonds [mol]))

(defrecord Molecule [_graph _name]
  PHasAtoms
  (atoms [mol] (nodes _graph))

  PMolecule
  (bonds [mol] (edges _graph))

  PMass
  (mass [mol] (reduce + (map mass (atoms mol))))
  (exact-mass [mol] (reduce + (map exact-mass (atoms mol))))

  clojure.lang.Named
  (getName [mol] _name))

(defn make-atom
  ([element name]
     (Atom. name (element/get-element element) nil nil 0 nil nil))
  ([element name isotope chirality charge hybridization explicit-hydrogen-count]
     (Atom. name (element/get-element element)
            isotope chirality charge hybridization explicit-hydrogen-count)))

(defn make-bond [atom1 atom2 & {:keys [type order direction]}]
  (Bond. [atom1 atom2] type order direction))

(defn make-molecule
  ([]
     (Molecule. (make-graph) nil))
  ;;; this is broken!!!
  ([atoms atom-pairs-vec]
     (reduce
      (fn [mol [atom1 atom2]]
        (assoc mol :_graph (add-edge (:_graph mol) atom1 atom2 (make-bond atom1 atom2))))
      (Molecule. (make-graph atoms) nil)
      atom-pairs-vec))
  ([atoms atom-pairs-vec name]
     (assoc (make-molecule atoms atom-pairs-vec) :name name)))

(defn add-atom [mol atom]
  (assoc mol :_graph (add-node (:_graph mol) atom)))

(defn get-atom [mol name]
  (get (atoms mol) name))

(defn add-bond
  ([mol bond]
     (assoc mol :_graph (add-edge (:_graph mol) bond)))
  ([mol atom1 atom2 & args]
     (assoc mol :_graph (add-edge (:_graph mol)
                                  (apply make-bond atom1 atom2 args)))))

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
