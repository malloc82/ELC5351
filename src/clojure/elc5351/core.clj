(ns elc5351.core
  (:use clojure.core
        elc5351.plot_patch)
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.operators :as op]
            [clojure.core.matrix.linear    :refer [norm]]
            [incanter.core   :as incanter]
            [incanter.charts :as charts]
            [incanter.stats  :as stats])
  (:import [mikera.matrixx Matrix Matrixx]
           [mikera.matrixx.impl SparseColumnMatrix SparseRowMatrix]
           [mikera.vectorz AVector Vector Vectorz]
           [mikera.vectorz.impl StridedVector SparseIndexedVector ZeroVector ArraySubVector AxisVector]
           [mikera.indexz Index]
           [java.awt GraphicsEnvironment]
           [java.util Arrays ArrayList Collection]))
(set! *warn-on-reflection* true)

(m/set-current-implementation :vectorz)


;; (def ^Matrix A (Matrix/createRandom 10 10))

(defn error-fn [^AVector a ^AVector b]
  (norm (m/sub a b)))

(defn art-error [^Matrix A x b]
  (let [Ax (m/inner-product A x)]
    (m/sub! Ax b)
    (norm Ax)))

(defn matrix-error [^Matrix A ^Matrix A_inv]
  (let [[rows cols]    (m/shape A)
        I      ^Matrix (m/identity-matrix rows)
        AA_inv ^Matrix (m/mmul A A_inv)
        error  ^AVector (Vector/createLength rows)]
    (dotimes [i rows]
      (m/mset! error i (norm (m/sub (m/get-column I i)
                                 (m/get-column AA_inv i)))))
    error))

(defn art-update [^Matrix A ^AVector b ^AVector x ^Double lambda & {:keys [index-order]}]
  (loop [indices (or index-order (range (first (m/shape A))))]
    (if indices
      (let [i ^Integer (int (first indices))
            A_i ^AVector (m/get-row A i)]
        (let [mag ^Double (.get (.innerProduct A_i A_i))
              diff ^Double (- (.unsafeGet b i)
                              (.get (.innerProduct A_i x)))]
          (.addMultiple x A_i ^Double (double (* lambda (/ diff mag))))
          (recur (next indices))))
      x)))

(defn art [^Matrix A ^AVector b & {:keys [iterations lambda verbose win series]
                                :or {iterations 10
                                     lambda 0.01}}]
  (let [[rows cols] (m/shape A)]
   (loop [i (int 0)
          x ^AVector (Vectorz/createUniformRandomVector cols)]
     (if (< i iterations)
       (let [new-x (art-update A b x lambda)
             error (art-error A new-x b)]
         (when verbose
           (println (format "test: iteration = %3d, error = %.6f" i error))
           (flush))
         (when win
           (update-series win series i error))
         (recur (unchecked-inc-int i) new-x))
       (do (when verbose
             (println "done.")
             (m/pm x))
           x)))))

(defn ^Matrix find-inverse [^Matrix A & {:keys [iterations lambda verbose plot-errors]
                                         :or {iterations 100
                                              lambda 0.5}}]
  (let [[rows cols] (m/shape A)
        I ^Matrix (m/identity-matrix rows)
        A_inverse ^Matrix (Matrix/create cols rows)
        win (when plot-errors
              (let [ch (charts/xy-plot [] []
                                       :series-label "default"
                                       :x-label "Iteration"
                                       :y-label "Error"
                                       :title (format "Matrix (%d X %d) Inverse Error Plots [lambda = %.3f, iterations = %d]"
                                                      rows cols lambda iterations)
                                       :legend true)]
                (remove-series ch "default")
                (incanter/view ch)
                ch))]
    (dotimes [i cols]
      (let [series-label (format "column %d" i)]
        (print ". ")
        (flush)
        (when win
          (charts/add-lines win [] [] :series-label series-label))
        (m/assign! (m/get-column A_inverse i)
                   (art A (m/get-column I i)
                        :iterations iterations
                        :lambda lambda
                        :verbose verbose
                        :win win
                        :series series-label))))
    (print "\n") (flush)
    A_inverse))

