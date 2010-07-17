(ns ugraph)

(defprotocol Graph
  (nodes [this])
  (edges [this] [this node])
  (edge? [this node1 node2])
  (neighbors [this node]))

(defprotocol Edge
  (nodes [this])
  (contains [this y])
  (other-node [this y]))

(defrecord UndirectedEdge [node1 node2]
  Edge
  (nodes [this] (vector node1 node2))
  (contains [this node] (or (= node node1) (= node node2)))
  (other-node [this node] (cond
                           (= node node1) node2
                           (= node node2) node1)))

(defrecord UndirectedGraph [nodes edges]
  Graph
  (edges [g]
         (distinct (apply concat (map vals (vals (:edges g))))))
  (edges [g node]
         (vals (get (:edges g) node)))
  (edge? [g n1 n2]
         (some #(contains % n2) (vals (get (:edges g) n1))))
  (neighbors [g node]
             (map #(other-node % node) (vals (get (:edges g) node)))))

(defn make-ugraph
  ([] (UndirectedGraph. #{} {}))
  ([nodes] (UndirectedGraph. nodes {}))
  ([nodes edges] (UndirectedGraph. nodes edges)))

(defn add-node [g n]
  (make-ugraph (conj (:nodes g) n) (:edges g)))

(letfn [(add-1-edge [e n1 n2 obj]
                    (assoc e n1 (assoc (or (get e n1) {}) n2 obj)))]
  (defn add-edge
    ([g n1 n2]
       (add-edge g n1 n2 (UndirectedEdge. n1 n2)))
    ([g n1 n2 obj]
       (let [e (:edges g)]
         (if (some #(contains % n2) (edges g n1))
           g
           (make-ugraph (:nodes g)
                        (add-1-edge
                         (add-1-edge e n2 n1 obj)
                         n1 n2 obj)))))))


(comment (defn bfs
           ([g start ] (bfs g start #{}))
            ([g start visited]
               (letfn
                   [(bfs-visit
                     [nodes-to-visit]
                     (map (fn [node])
                          nodes-to-visit))]
                 (bfs-visit (list start))))))