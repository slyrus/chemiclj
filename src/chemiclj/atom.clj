;;; file: chemiclj/atom.clj
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


(ns chemiclj.atom
  (:use [chemiclj.mass])
  (:require [chemiclj.element :as element]))

(defprotocol PAtomContainer
  (atoms [obj]))

;;;
;;; Atom record and related functions
(defrecord Atom [_name element isotope chirality charge
                 hybridization aromatic explicit-hydrogen-count]
  chemiclj.mass.PMass
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

