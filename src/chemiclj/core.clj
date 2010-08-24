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


;;;
;;; CORE NAMESPACE STARTS HERE!!!!
;;;
(ns chemiclj.core
  (:use [chemiclj element atom bond molecule])
  (:require [chemiclj.protocol :as protocol]
            [clojure.contrib.def :as d]
            [shortcut.graph :as g]))

(d/defalias neighbors g/neighbors)

(d/defalias isotopes protocol/isotopes)

(d/defalias mass protocol/mass)
(d/defalias exact-mass protocol/exact-mass)

(d/defalias atoms protocol/atoms)
(d/defalias get-atom protocol/get-atom)

(d/defalias add-atom protocol/add-atom)
(d/defalias remove-atom protocol/remove-atom)
(d/defalias bonds protocol/bonds)
(d/defalias bond? protocol/bond?)
(d/defalias add-bond protocol/add-bond)
(d/defalias remove-bond protocol/remove-bond)
(d/defalias configurations protocol/configurations)
(d/defalias add-configuration protocol/add-configuration)

(defn names [seq]
  (map name seq))

(defn molecular-formula [mol]
  (apply str
         (map #(str (:id (first %)) (second %))
              (sort-by #(:id (first %)) (count-elements mol)))))

(defn get-atoms-of-element [mol elmnt]
  (let [elmnt (get-element elmnt)]
    (filter #(= (:element %) elmnt) (atoms mol))))

(defn remove-atoms-of-element [mol element]
  (reduce (fn [mol atom]
            (remove-atom mol atom))
          mol
          (get-atoms-of-element mol element)))

