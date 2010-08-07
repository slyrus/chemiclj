;;; file: chemiclj/element.clj
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

(ns chemiclj.element
  (:require [clojure.string :as string]
            [clojure.xml :as xml]
            [clojure.zip :as zip]
            [clojure.contrib.zip-filter.xml :as zf]))

(defprotocol HasIsotopes
  (isotopes [obj]))

(defrecord Isotope [element number id exact-mass relative-abundance])

(defn element-abundnant-isotopes [element]
  (map second (reverse (sort-by #(:relative-abundance (val %))
                                (isotopes element)))))

(defn get-scalar-text [node dict-ref]
  (zf/xml1->
   (zip/xml-zip node)
   :scalar [(zf/attr= :dictRef dict-ref)]
   zf/text))

(defn float-if [number?]
  (when number?
    (Float/parseFloat number?)))

(defn int-if [number?]
  (when number?
    (Integer/parseInt number?)))

(def isotope-list
     (zf/xml->
      (zip/xml-zip
       (xml/parse "data/isotopes.xml"))
      :isotopeList
      zip/right
      zip/node
      (fn [isotope-list]
        (let [{:keys [id]} (:attrs isotope-list)]
          (zf/xml->
           (zip/xml-zip isotope-list)
           :isotope
           zip/node
           (fn [isotope]
             (let [{:keys [id number elementType]} (:attrs isotope)
                   relative-abundance (get-scalar-text isotope "bo:relativeAbundance")
                   exact-mass (get-scalar-text isotope "bo:exactMass")]
               (Isotope. elementType (int-if number) id
                         (float-if exact-mass)
                         (float-if relative-abundance)))))))))

(defn convert-alist
  [alist & [key first]]
  (reduce (fn [v x] (assoc v (key x) x))
          (vec (repeat (apply max (map key alist)) nil))
          alist))

(defrecord Element [atomic-number id name group period mass
                    electronegativity max-bond-order]
  HasIsotopes
  (isotopes [element]
            (reduce (fn [v x] (assoc v (:number x) x))
                    (hash-map)
                    (filter (fn [x] (= (:element x) id))
                            isotope-list)))
  clojure.lang.Named
  (getName [atm] name))

(def elements-list (zf/xml->
                    (zip/xml-zip
                     (xml/parse "data/elementdata.xml"))
                    zip/down
                    zip/rights
                    (fn [x]
                      (let [{:keys [atomicnumber id name group period]} (:attrs x)
                            mass (zf/xml1-> (zip/xml-zip x) :mass zf/text)
                            maxbondorder
                            (zf/xml1-> (zip/xml-zip x)
                                       :maxbondorder zf/text)
                            electronegativity
                            (zf/xml1-> (zip/xml-zip x)
                                       :electronegativity zf/text)]
                        (Element. (Integer/parseInt atomicnumber)
                                  id
                                  name
                                  (Integer/parseInt group)
                                  (Integer/parseInt period)
                                  (Float/parseFloat mass)
                                  (when electronegativity
                                    (Float/parseFloat electronegativity))
                                  (when maxbondorder
                                    (Integer/parseInt maxbondorder))
                                  )))))

(def elements-vec (convert-alist elements-list :atomic-number))

(def element-id-hash (reduce (fn [v x] (assoc v (:id x) x))
                             (hash-map)
                             elements-vec))

(def element-name-hash (reduce (fn [v x] (assoc v (:name x) x))
                               (hash-map)
                               elements-vec))

(defmulti get-element class)

(defmethod get-element Number [id] (nth elements-vec id))
(defmethod get-element String [id]
  (or
   (get element-id-hash (str (string/upper-case (nth id 0))
                             (string/lower-case (subs id 1 (count id)))))
   (get element-name-hash (string/lower-case id))))
(defmethod get-element clojure.lang.Keyword [id]
  (get-element (name id)))

(def *element-normal-valences*
     (reduce (fn [v x] (assoc v (get-element (first x)) (rest x)))
             (hash-map)
             '(("B" 3)
               ("C" 4)
               ("N" 3 5)
               ("O" 2)
               ("P" 3 5)
               ("S" 2 4 6)
               ("F" 1)
               ("Cl" 1)
               ("Br" 1)
               ("I" 1))))

(defmulti get-normal-valences class)

(defmethod get-normal-valences Number [id] (get *element-normal-valences* (get-element id)))
(defmethod get-normal-valences String [id] (get *element-normal-valences* (get-element id)))
(defmethod get-normal-valences Element [element] (get *element-normal-valences* element))

