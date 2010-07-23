
(ns chemiclj.core
  (:use clojure.java.io)
  (:use shortcut.graph)

  (:require [clojure.string :as string]
            [clojure.xml :as xml]
            [clojure.zip :as zip]
            [clojure.contrib.zip-filter.xml :as zf])

  (:import shortcut.graph.Graph))

(defrecord Element [atomic-number id name group period mass
                    electronegativity max-bond-order isotopes])

(defrecord Isotope [element number id exact-mass relative-abundance])

(defrecord Atom [name element charge isotope hybridization])

(defprotocol BondOps
  (atom-count [bond]))

(defrecord Bond [nodes type order direction]
  NodeSet
  (nodes [bond] nodes)
  (node? [bond node] (some #{node} nodes))
  (neighbors [bond node] (remove #{node} nodes))

  clojure.lang.Indexed
  (nth [bond i] (nth nodes i))

  BondOps
  (atom-count [bond] (count nodes)))

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
       (xml/parse "/Users/sly/projects/chemiclj/data/isotopes.xml"))
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

(def elements-list (zf/xml->
                    (zip/xml-zip
                     (xml/parse "/Users/sly/projects/chemiclj/data/elementdata.xml"))
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
                                  (reduce (fn [v x] (assoc v (:number x) x))
                                          (hash-map)
                                          (filter (fn [x] (= (:element x) id))
                                                  isotope-list)))))))

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

(defprotocol Molecule
  (mass [mol])
  (exact-mass [mol]))

(defn element-abundnant-isotopes [element]
  (map second (reverse (sort-by #(:relative-abundance (val %))
                                (:isotopes element)))))

(defn atom-exact-mass [atom]
  (if (:isotope atom)
    (:exact-mass (:isotope atom))
    (:exact-mass (first (element-abundnant-isotopes (:element atom))))))

(defn mass [mol]
  (reduce + (map #(-> % :element :mass) (nodes mol))))

(defn exact-mass [mol]
  (reduce + (map atom-exact-mass (nodes mol))))

(defn make-atom
  ([element name]
     (Atom. name (get-element element)
            nil
            nil
            nil)))

(defn make-bond [atom1 atom2 & {:keys [type order direction],
                                :or {type :single order 1}}]
  (Bond. [atom1 atom2] type order direction))

(defn make-molecule [atoms bond-vec]
  (make-graph atoms (vec (map (fn [[a b]] (make-bond a b)) bond-vec))))
