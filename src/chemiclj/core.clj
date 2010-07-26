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
  (:use chemiclj.element)

  (:use shortcut.graph)

  (:require [clojure.string :as string]))

(defprotocol PMass
  (mass [mol])
  (exact-mass [mol]))

(defrecord Atom [name element charge isotope hybridization]
  PMass
  (mass [atm] (-> atm :element :mass))
  (exact-mass [atm]
              (if (:isotope atm)
                (:exact-mass (:isotope atm))
                (:exact-mass (first (element-abundnant-isotopes (:element atm)))))))

(defrecord Bond [nodes type order direction]
  NodeSet
  (nodes [bond] nodes)
  (node? [bond node] (some #{node} nodes))
  (neighbors [bond node] (remove #{node} nodes))

  clojure.lang.Indexed
  (nth [bond i] (nth nodes i)))

(defprotocol PMolecule
  (atoms [mol])
  (bonds [mol]))

(defrecord Molecule [_graph _name]

  PMolecule
  (atoms [mol] (nodes _graph))
  (bonds [mol] (edges _graph))

  PMass
  (mass [mol] (reduce + (map mass (atoms mol))))
  (exact-mass [mol] (reduce + (map exact-mass (atoms mol)))))

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
