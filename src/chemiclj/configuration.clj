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

(ns chemiclj.configuration
  (:require [clojure.contrib [except :as except]]))

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


