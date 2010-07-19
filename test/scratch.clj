
(ns chemiclj-test
  (:use chemiclj.core)
  (:use shortcut.graph)

  (:import chemiclj.core.Atom)
  (:import shortcut.graph.Graph))

(defn make-atom [element name]
  (Atom. name (get-element element) nil nil nil))

(def c1 (make-atom :c "C1"))
(def h1 (make-atom :h "H1"))
(def h2 (make-atom :h "H2"))
(def h3 (make-atom :h "H3"))
(def c2 (make-atom :c "C2"))
(def h4 (make-atom :h "H4"))
(def h5 (make-atom :h "H5"))
(def h6 (make-atom :h "H6"))

(def ethane (add-edges (Graph. #{c1 h1 h2 h3 c2 h4 h5 h6} {})
                       [[c1 h1]
                        [c1 h2]
                        [c1 h3]
                        [c1 c2]
                        [c2 h4]
                        [c2 h5]
                        [c2 h6]]))
