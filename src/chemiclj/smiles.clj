
(ns chemiclj.smiles
  (:use [chemiclj.core])
  (:require [chemiclj.element :as element]
            [shortcut.graph :as graph])
  
  (:import chemiclj.core.Atom))

(defn collect-while [pred s]
  "similar to take-while, except that predicate takes no arguments and
the returned sequence isn't lazy."
  (when (pred)
    (cons (first s)
          (collect-while pred (rest s)))))

(def *aromatic-atoms* ["c" "n" "o" "p" "s" "as" "se"])
(def *organic-atoms* ["B" "C" "N" "O" "P" "S" "F" "Cl" "Br" "I"])

(defn parse-smiles [string & [name add-implicit-hydrogens]]
  (let [pos (atom 0)
        counts (atom {})]
    (letfn
        [(peek-char [] (and (< @pos (count string))
                            (nth string @pos)))
         (read-char [] (let [c (nth string @pos)]
                         (swap! pos inc)
                         c))
         (skip [n] (swap! pos + n))
         (add-atom [mol element last]
                   (let [ind (inc (or (get @counts element) 0))]
                     (swap! counts assoc element ind)
                     (let [atom (make-atom element (str (:id (element/get-element element)) ind))]
                       [(if last
                          (add-bond (chemiclj.core/add-atom mol atom) last atom)
                          (chemiclj.core/add-atom mol atom))
                        atom])))
         (read-number []
                      (let [num (collect-while
                                 #(let [c (peek-char)]
                                    (and c (not (neg? (Character/digit c 10)))))
                                 (repeatedly #(Character/digit (read-char) 10)))]
                        (when (seq num)
                          (Integer/parseInt (apply str num)))))
         (read-bracket-expression [mol last])
         (read-branch [mol last])
         (read-smiles-token [mol last]
                            (when (peek-char)
                              (let [c (read-char)]
                                (cond
                                 (= c \[) (read-bracket-expression mol last)
                                 (= c \() (read-branch mol last)

                                 (= c \-) (list :bond :single)
                                 (= c \=) (list :bond :double)
                                 (= c \#) (list :bond :triple)
                                 (= c \:) (list :bond :aromatic)
                                 (= c \/) (list :bond :up)
                                 (= c \\) (list :bond :down)

                                 (or (and c (not (neg? (Character/digit c 10))))
                                     (= c \%))
                                 (let [number (read-number)]
                                   (print number))
                               
                                 (= c \B)
                                 (if (= (peek-char) \r)
                                   (do (read-char)
                                       (add-atom mol "Br" last))
                                   (add-atom mol "B" last))
                               
                                 (= c \C)
                                 (if (= (peek-char) \l)
                                   (do (read-char)
                                       (add-atom mol "Cl" last))
                                   (add-atom mol "C" last))

                                 (= c \N) (add-atom mol "N" last)
                                 (= c \O) (add-atom mol "O" last)
                                 (= c \P) (add-atom mol "P" last)
                                 (= c \S) (add-atom mol "S" last)
                                 (= c \F) (add-atom mol "F" last)
                                 (= c \I) (add-atom mol "I" last)))))]
      (loop [mol (apply make-molecule (when name (list name)))
             last nil]
        (let [v (read-smiles-token mol last)]
          (if v
            (let [[mol new-last] v]
              (recur mol (or new-last last)))
            mol))))))


(defn foo [string]
  (let [pos (atom 0)]
    (letfn
        [(peek-char [] (and (< @pos (count string))
                            (nth string @pos)))
         (read-char [] (let [c (nth string @pos)]
                         (swap! pos inc)
                         c))
         (skip [n] (swap! pos + n))]
      [(read-char)
       (peek-char)])))