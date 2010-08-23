;;; file: chemiclj/molecule.clj
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

(ns chemiclj.molecule
  (:use [chemiclj [core atom bond]])
  (:require [shortcut.graph :as g]
            [clojure.contrib.def :as d]
            [clojure.contrib [except :as except]]))

;;;
;;; Molecule record and related functions
(defrecord Molecule [_graph _name _configurations]
  chemiclj.core.PAtomContainer
  (atoms [mol] (g/nodes _graph))

  chemiclj.core.PMolecule
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

