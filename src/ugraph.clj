;;; file: ugraph.clj
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

(ns ugraph)

(defprotocol NodeSet
  (nodes [graph])
  (node? [graph node])
  (add-node [graph node])
  (neighbors [graph node]))

(defprotocol EdgeSet
  (edges [graph]
         [graph node])
  (edge? [graph node1 node2])
  (add-edge [graph node1 node2]
            [graph node1 node2 object]))

(defprotocol ArcInGraph
  (start [edge])
  (end [edge]))

(defrecord Edge [node1 node2]
  NodeSet
  (nodes [edge] (vector node1 node2))
  (node? [edge node] (or (= node node1) (= node node2)))
  (neighbors [edge node] (list (cond
                                (= node node1) node2
                                (= node node2) node1))))

(defrecord Arc [start-node end-node]
  NodeSet
  (nodes [edge] (vector start-node end-node))
  (node? [edge node] (or (= node start-node) (= node end-node)))
  (neighbors [edge node] (list (cond
                                (= node start-node) end-node
                                (= node end-node) start-node)))
  ArcInGraph
  (start [edge] start-node)
  (end [edge] end-node))

;;; forward declaration for make-ugraph which we can't defn until we
;;; defrecord Graph
(declare make-ugraph)

(defrecord Graph [node-set edge-map]
  NodeSet
  (nodes [g] (:node-set g))
  (node? [g node] (get (nodes g) node))
  (add-node [g n]
            (make-ugraph (conj (:node-set g) n) (:edge-map g)))
  (neighbors [g node]
             (map #(first (neighbors % node)) (vals (get (:edge-map g) node))))
  
  EdgeSet
  (edges [g]
         (distinct (apply concat (map vals (vals (:edge-map g))))))
  (edges [g node]
         (vals (get (:edge-map g) node)))
  
  (edge? [g n1 n2]
         (some #(when (node? % n2) %)
               (vals (get (:edge-map g) n1))))
  (add-edge [g n1 n2]
            (add-edge g n1 n2 (Edge. n1 n2)))
  (add-edge [g n1 n2 obj]
            (letfn [(add-1-edge [e n1 n2 obj]
                                (assoc e n1 (assoc (or (get e n1) {}) n2 obj)))]
              (if (some #(node? % n2) (edges g n1))
                g
                (make-ugraph (:node-set g)
                             (add-1-edge
                              (add-1-edge (:edge-map g) n2 n1 obj)
                              n1 n2 obj))))))

(defn make-ugraph
  ([] (Graph. #{} {}))
  ([nodes] (Graph. nodes {}))
  ([nodes edges] (Graph. nodes edges)))

(defn add-edges [g edge-vec]
  (reduce (fn [g [n1 n2]] (add-edge g n1 n2)) g edge-vec))

(defn add-nodes [g & node-vec]
  (reduce (fn [g node] (add-node g node)) g node-vec))

(defn breadth-first-traversal
  ([g start]
     (when (node? g start)
       (breadth-first-traversal g [start] #{})))
  ([g queue visited]
     (lazy-seq
      (if (seq queue)
        (let [node (first queue)
              next (remove visited (neighbors g node))]
          (cons node
                (breadth-first-traversal g (into (subvec queue 1) next)
                                         (into (conj visited node) next))))))))

;;; often it's nice to not just do a search, but keep a trail of the
;;; path of how one got to a particular node. we're also going to want
;;; to know the distance from the start, so we can just take the
;;; length of the path to get that.
(defn breadth-first-traversal-with-path
  ([g start]
     (when (node? g start)
       (breadth-first-traversal-with-path g [[start]] #{})))
  ([g queue visited]
     (lazy-seq
      (if (seq queue)
        (let [node (last (first queue))
              next (remove visited (neighbors g node))]
          (cons (first queue)
                (breadth-first-traversal-with-path
                  g
                  (into (subvec queue 1)
                        (vec (map #(conj (first queue) %) (vec next))))
                  (into (conj visited node) next))))))))

(defn depth-first-traversal
  ([g start]
     (when (node? g start)
       (depth-first-traversal g (list start) #{})))
  ([g queue visited]
     (lazy-seq
      (if (seq queue)
        (let [node (first queue)
              next (remove visited (neighbors g node))]
          (if-not (visited node)
            (cons node
                  (depth-first-traversal g (into (rest queue) next)
                                         (conj visited node)))
            (depth-first-traversal g (into (rest queue) next) visited)))))))

(defn depth-first-traversal-with-path
  ([g start]
     (when (node? g start)
       (depth-first-traversal-with-path g (list [start]) #{})))
  ([g queue visited]
     (lazy-seq
      (if (seq queue)
        (let [node (last (first queue))
              next (remove visited (neighbors g node))]
          (if-not (visited node)
            (cons (first queue)
                  (depth-first-traversal-with-path g
                    (into (rest queue)
                          (vec (map #(conj (first queue) %) (vec next))))
                    (conj visited node)))
            (depth-first-traversal-with-path g
                        (into (rest queue)
                              (vec (map #(conj (first queue) %) (vec next))))
                        visited)))))))


;;; scratch

(def q (make-ugraph #{1 2 3 4} {}))
(def q2 (add-edge (add-edge q 3 1) 1 4))
(def q3 (add-edges q [[1 2] [1 3] [3 4]]))
(def q4 (add-edges (add-nodes q3 5 6) [[4 5] [1 5] [4 6]]))
(def q5 (add-edges (make-ugraph #{1 2 3 4 5 6})
                   [[1 2] [2 3] [3 4] [4 5] [5 6] [6 1]]))

(def q6 (reduce #(add-node %1 %2) (make-ugraph) (range 1000)))
(def q7 (reduce (fn [g [n1 n2]] (add-edge g n1 n2))
                q6
                (take 10000 (repeatedly #(vector (rand-int 1000)
                                                 (rand-int 1000))))))
(neighbors q7 1)
(take-while (complement #{2}) (breadth-first-traversal q7 1))

(def q8 (reduce #(add-node %1 (str "node" %2)) (make-ugraph) (range 1000)))


(def q9 (reduce (fn [g [n1 n2]] (add-edge g
                                          (str "node" n1)
                                          (str "node" n2)))
                q8
                (take 10000 (repeatedly #(vector (rand-int 1000)
                                                 (rand-int 1000))))))
