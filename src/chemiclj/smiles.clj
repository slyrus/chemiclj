
(ns chemiclj.smiles
  (:use [chemiclj.core])
  (:require [chemiclj.element :as element]
            [shortcut.graph :as graph]
            [clojure.contrib.error-kit :as err])
  
  (:import chemiclj.core.Atom))


(err/deferror *smiles-parser-error* [] [str]
  {:msg (str "Smiles Parser Error: " str)
   :unhandled (err/throw-msg IllegalArgumentException)})

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
         (add-atom [mol element state]
                   (let [ind (inc (or (get @counts element) 0))]
                     (swap! counts assoc element ind)
                     (let [atom
                           (make-atom element
                                      (str (:id (element/get-element element)) ind))]
                       [(if state
                          (add-bond (chemiclj.core/add-atom mol atom) state atom)
                          (chemiclj.core/add-atom mol atom))
                        atom])))
         (read-number []
                      (let [num (collect-while
                                 #(let [c (peek-char)]
                                    (and c (not (neg? (Character/digit c 10)))))
                                 (repeatedly #(Character/digit (read-char) 10)))]
                        (when (seq num)
                          (Integer/parseInt (apply str num)))))
         (read-bracket-expression [mol state])
         (read-branch [mol state]
                      (let [mol (read-smiles-tokens mol state)]
                        [mol state]))
         (read-smiles-token [mol state]
                            (when (peek-char)
                              (let [c (read-char)]
                                (cond
                                 (= c \[) (read-bracket-expression mol state)

                                 (= c \() (read-branch mol state)
                                 (= c \)) nil

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
                                       (add-atom mol "Br" state))
                                   (add-atom mol "B" state))
                               
                                 (= c \C)
                                 (if (= (peek-char) \l)
                                   (do (read-char)
                                       (add-atom mol "Cl" state))
                                   (add-atom mol "C" state))

                                 (= c \N) (add-atom mol "N" state)
                                 (= c \O) (add-atom mol "O" state)
                                 (= c \P) (add-atom mol "P" state)
                                 (= c \S) (add-atom mol "S" state)
                                 (= c \F) (add-atom mol "F" state)
                                 (= c \I) (add-atom mol "I" state)))))
         (read-smiles-tokens [mol state]
                             (let [v (read-smiles-token mol state)]
                               (if v
                                 (let [[mol new-state] v]
                                   (recur mol (or new-state state)))
                                 mol)))]
      (loop [mol (apply make-molecule (when name (list name)))
             state nil]
        (read-smiles-tokens mol state)))))


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