
(ns chemiclj
  (:use clojure.java.io)
  (:use clojure.contrib.graph)

  (:require [clojure.xml :as xml]
            [clojure.zip :as zip]
            [clojure.contrib.zip-filter.xml :as zf]))

(defrecord Element [atomic-number id name group period mass
                    electronegativity max-bond-order])

(def elements-zip
     (zip/xml-zip (xml/parse "/Users/sly/projects/chemiclj/data/elementdata.xml")))

(def elements-list
     (zf/xml-> elements-zip zip/down zip/rights
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
                               (Integer/parseInt maxbondorder)))))))

(defn convert-alist
  [alist & [key first]]
  (reduce (fn [v x] (assoc v (key x) x))
          (vec (repeat (apply max (map key alist)) nil))
          alist))

(def elements-vec (convert-alist elements-list :atomic-number))

(defrecord Isotope [number exact-mass relative-abundance])

(defrecord Atom [name element charge isotope hybridization])

