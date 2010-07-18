
(ns ugraph)

(defprotocol NodeSet
  (nodes [graph]))

(defprotocol Graph
  (add-node [graph node])
  (add-edge [graph node1 node2]
            [graph node1 node2 object])
  (edges [graph]
         [graph node])
  (edge? [graph node1 node2])
  (neighbors [graph node]))

(defprotocol Edge
  (contains [edge y])
  (other-node [edge y]))

(defrecord UndirectedEdge [node1 node2]
  NodeSet
  (nodes [edge] (vector node1 node2))
  Edge
  (contains [edge node] (or (= node node1) (= node node2)))
  (other-node [edge node] (cond
                           (= node node1) node2
                           (= node node2) node1)))

;;; forward declaration for make-ugraph which we can't defn until we
;;; defrecord UndirectedGraph
(declare make-ugraph)

(defrecord UndirectedGraph [node-set edge-map]
  NodeSet
  (nodes [g] (:node-set g))
  Graph
  (add-node [g n]
            (make-ugraph (conj (:node-set g) n) (:edge-map g)))
  (edges [g]
         (distinct (apply concat (map vals (vals (:edge-map g))))))
  (edges [g node]
         (vals (get (:edge-map g) node)))
  (edge? [g n1 n2]
         (some #(when (contains % n2) %)
               (vals (get (:edge-map g) n1))))
  (add-edge [g n1 n2]
            (add-edge g n1 n2 (UndirectedEdge. n1 n2)))
  (add-edge [g n1 n2 obj]
            (letfn [(add-1-edge [e n1 n2 obj]
                                (assoc e n1 (assoc (or (get e n1) {}) n2 obj)))]
              (if (some #(contains % n2) (edges g n1))
                g
                (make-ugraph (:node-set g)
                             (add-1-edge
                              (add-1-edge (:edge-map g) n2 n1 obj)
                              n1 n2 obj)))))
  (neighbors [g node]
             (map #(other-node % node) (vals (get (:edge-map g) node)))))

(defn make-ugraph
  ([] (UndirectedGraph. #{} {}))
  ([nodes] (UndirectedGraph. nodes {}))
  ([nodes edges] (UndirectedGraph. nodes edges)))

(defn add-edges [g edge-vec]
  (reduce (fn [g [n1 n2]] (add-edge g n1 n2)) g edge-vec))

(defn add-nodes [g & node-vec]
  (reduce (fn [g node] (add-node g node)) g node-vec))

(defn breadth-first-traversal
  ([g start]
     (breadth-first-traversal g [start] #{}))
  ([g queue visited]
     (lazy-seq
      (if (seq queue)
        (let [node (first queue)
              next (remove visited (neighbors g node))]
          (cons node (breadth-first-traversal g (into (subvec queue 1) next)
                                              (into (conj visited node) next))))))))

(defn depth-first-traversal
  ([g start]
     (depth-first-traversal g (list start) [[] #{}]))
  ([g queue [path visited]]
     (lazy-seq
      (if (seq queue)
        (let [node (first queue)
              next (remove visited (neighbors g node))]
          (depth-first-traversal g (into (rest queue) next)
                                 [(conj path node)
                                  (into (conj visited node) next)]))
        path))))


;;; scratch

(def q (make-ugraph #{1 2 3 4} {}))
(def q2 (add-edge (add-edge q 3 1) 1 4))
(def q3 (add-edges q [[1 2] [1 3] [3 4]]))
(def q4 (add-edges (add-nodes q3 5 6) [[4 5] [1 5] [4 6]]))
(def q5 (add-edges (make-ugraph #{1 2 3 4 5 6})
                   [[1 2] [2 3] [3 4] [4 5] [5 6] [6 1]]))
