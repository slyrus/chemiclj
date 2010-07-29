
(ns chemiclj.smiles
  (:use [chemiclj.core])
  (:require [shortcut.graph :as graph])
  
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
  (let [pos (atom 0)]
    (letfn
        [(peek-char [] (and (< @pos (count string))
                            (nth string @pos)))
         (read-char [] (let [c (nth string @pos)]
                         (swap! pos inc)
                         c))
         (skip [n] (swap! pos + n))
         (add-atom [mol element counter]
                   (let [atm (make-atom element counter)]
                     (chemiclj.core/add-atom mol atm)))
         (read-number []
                      (let [num (collect-while
                                 #(let [c (peek-char)]
                                    (and c (not (neg? (Character/digit c 10)))))
                                 (repeatedly #(Character/digit (read-char) 10)))]
                        (when (seq num)
                          (Integer/parseInt (apply str num)))))
         (read-bracket-expression [])
         (read-branch [])
         (read-atom [])
         (read-smiles-token [mol last counter]
                            (when (peek-char)
                              (let [c (read-char)]
                                (cond
                                 (= c \[) (read-bracket-expression)
                                 (= c \() (read-branch)
                                 (= c \)) nil
                                 (= c \-) (list :bond :single)
                                 (= c \=) (list :bond :double)
                                 (= c \#) (list :bond :triple)
                                 (= c \:) (list :bond :aromatic)
                                 (= c \/) (list :bond :up)
                                 (= c \\) (list :bond :down)
                                 (= c \.) (list :disconnected)

                                 (or (and c (not (neg? (Character/digit c 10))))
                                     (= c \%))
                                 (let [number (read-number)]
                                   (print number))
                               
                                 (= c \B)
                                 (if (= (peek-char) \r)
                                   (do (read-char)
                                       (list :atom :bromine))
                                   (list :atom :boron))
                               
                                 (= c \C)
                                 (if (= (peek-char) \l)
                                   (do (read-char)
                                       (list :atom :chlorine))
                                   (list :atom :carbon))

                                 (= c \N) [(add-atom mol "N" counter)]
                                 (= c \O) [(add-atom mol "O" counter)]
                                 (= c \P) (list :atom :phosphorus)
                                 (= c \S) (list :atom :sulfur)
                                 (= c \F) (list :atom :fluorine)
                                 (= c \I) (list :atom :iodine)))))]
      (loop [mol (apply make-molecule (when name (list name)))
             last nil
             counter 0]
        (let [v (read-smiles-token mol last counter)]
          (if v
            (let [[mol] v]
              (recur mol last (inc counter)))
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