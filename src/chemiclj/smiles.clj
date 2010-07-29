
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
         (add-atom [mol element state]
                   (let [ind (inc (or (get @counts element) 0))]
                     (swap! counts assoc element ind)
                     (let [atom
                           (make-atom element
                                      (str (:id (element/get-element element)) ind))]
                       [(let [last-atom (:last-atom state)]
                          (if last-atom
                            (let [bond-type (:next-bond state)]
                              (cond
                               (= bond-type :double)
                               (add-double-bond (chemiclj.core/add-atom mol atom) last-atom atom)

                               (= bond-type :triple)
                               (add-triple-bond (chemiclj.core/add-atom mol atom) last-atom atom)
                               
                               (= bond-type :aromatic)
                               (add-aromatic-bond (chemiclj.core/add-atom mol atom) last-atom atom)

                               true
                               (add-bond (chemiclj.core/add-atom mol atom) last-atom atom)))
                            (chemiclj.core/add-atom mol atom)))
                        (dissoc (assoc state :last-atom atom)
                                :next-bond)])))
         (read-number []
                      (let [num (collect-while
                                 #(let [c (peek-char)]
                                    (and c (not (neg? (Character/digit c 10)))))
                                 (repeatedly #(Character/digit (read-char) 10)))]
                        (when (seq num)
                          (Integer/parseInt (apply str num)))))
         (read-bracket-expression [mol state])
         (read-branch [mol state]
                      (let [[mol new-state]
                            (read-smiles-tokens mol state)]
                        ;; restore last atom before the branch as the
                        ;; last-atom in state so next atom is
                        ;; connected properly
                        [mol (assoc new-state :last-atom (:last-atom state))]))
         (read-smiles-token [mol state]
                            (when (peek-char)
                              (let [c (read-char)]
                                (cond
                                 (= c \[) (read-bracket-expression mol state)

                                 (= c \() (read-branch mol state)
                                 (= c \)) (list :end-branch state)

                                 (= c \-) [mol (assoc state :next-bond :single)]
                                 (= c \=) [mol (assoc state :next-bond :double)]
                                 (= c \#) [mol (assoc state :next-bond :triple)]
                                 (= c \:) [mol (assoc state :next-bond :aromatic)]
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
                               (cond
                                (nil? v) mol

                                (= (first v) :end-branch)
                                [mol (second v)]

                                true
                                (let [[mol new-state] v]
                                       (recur mol (or new-state state))))))]
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