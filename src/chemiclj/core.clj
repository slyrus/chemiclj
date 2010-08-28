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


(ns chemiclj.core)
  
(defprotocol PElement
  (isotopes [obj]))

(defprotocol PMass
  (mass [mol])
  (exact-mass [mol]))

(defprotocol PAtomContainer
  (atoms [obj])
  (get-atom [obj id]))

(defprotocol PMolecule
  (add-atom [mol atom])
  (remove-atom [mol atom])
  (bonds [mol] [mol atom])
  (bond? [mol atom1 atom2])
  (add-bond [mol bond] [mol atom1 atom2])
  (remove-bond [mol bond] [mol atom1 atom2])
  (configurations [mol])
  (add-configuration [mol configuration]))

(defmacro defn* [fn-name & rest]
  `(intern (if ~(namespace fn-name)
             (symbol ~(namespace fn-name))
             (ns-name *ns*))
           (symbol ~(name fn-name)) (fn ~(symbol (name fn-name)) ~@rest)))

(load "atom")
(load "molecule")
