(ns elc5351.utils
  (:use clojure.core
        elc5351.plot_patch)
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.operators :as op]
            [clojure.core.matrix.linear    :refer [norm]]
            [incanter.core   :as incanter]
            [incanter.charts :as charts]
            [incanter.stats  :as stats])
  (:import [mikera.matrixx Matrix Matrixx]
           [mikera.vectorz AVector Vector Vectorz]))

(defn ^Boolean sameShape? [^Matrix A ^Matrix B]
  (let [[rows_A cols_A] (m/shape A)
        [rows_B cols_B] (m/shape B)]
    (if (and (= rows_A rows_B) (= cols_A cols_B)) true false)))

(defn ^Matrix createComplexRandom
  ([^Integer rows]
   (createComplexRandom rows rows))
  ([^Integer rows ^Integer cols]
   (Matrix/createRandom  rows (bit-shift-left cols 1))))


(defn ^Matrix realToComplex
  ([^Matrix m]
   (let [[rows cols] (m/shape m)]
     (realToComplex m rows cols)))
  ([^Matrix m ^Integer rows ^Integer cols]
   (let [m_complex  ^Matrix  (Matrix/create rows (bit-shift-left cols 1))
         m_arr      ^doubles (.asDoubleArray m)
         m_comp_arr ^doubles (.asDoubleArray m_complex)]
     (dotimes [i (alength m_arr)]
       (aset m_comp_arr (bit-shift-left i 1) (aget m_arr i)))
     m_complex)))


(defn ^Matrix composeComplex [^Matrix amplitude ^Matrix phase & {:keys [target]}]
  (let [[rows_A cols_A] (m/shape amplitude)
        [rows_P cols_P] (m/shape phase)]
    (when-not (and (= rows_A rows_P) (= cols_A cols_P))
      (println (format "A : [%d, %d], B : [%d, %d]" rows_A cols_A rows_P cols_P))
      (throw (IndexOutOfBoundsException. "composeComplex: amplitude and phase have different dimension.")))
    (let [complex_mat ^Matrix  (Matrix/create rows_A (bit-shift-left cols_A 1))
          complex_arr ^doubles (.asDoubleArray complex_mat)
          amp_arr     ^doubles (.asDoubleArray amplitude)
          phase_arr   ^doubles (.asDoubleArray phase)]
      (dotimes [i (alength amp_arr)]
        (let [idx (bit-shift-left i 1)]
          (aset complex_arr idx (* (aget amp_arr i) (Math/cos (aget phase_arr i))))
          (aset complex_arr (unchecked-inc idx) (* (aget amp_arr i) (Math/sin (aget phase_arr i))))))
      complex_mat)))


(defn ^Matrix complexToReal [^Matrix m & {:keys [type]
                                          :or {type :mag}}]
  (let [[rows cols] (m/shape m)
        real_mat    ^Matrix  (Matrix/create rows (bit-shift-right cols 1))
        real_arr    ^doubles (.asDoubleArray real_mat)
        complex_arr ^doubles (.asDoubleArray m)]
    (case type
      :mag   (dotimes [i (alength real_arr)]
               (let [idx (bit-shift-left i 1)]
                 (aset real_arr i ^Double (norm [(aget complex_arr idx)
                                                 (aget complex_arr (unchecked-inc idx))]))))
      :phase (dotimes [i (alength real_arr)]
               (let [idx (bit-shift-left i 1)]
                 (aset real_arr i ^Double (Math/atan2 (aget complex_arr (unchecked-inc idx))
                                                      (aget complex_arr idx))))))
    real_mat))

(defn ^Double mat_norm [^Matrix m & {:keys [type] :or {type :inf}}]
  (let [[rows cols] (m/shape m)]
    (case type
      :inf (let [v ^AVector (Vector/createLength rows)
                 m_abs ^Matrix (m/abs m)]
             (dotimes [i rows]
               (.set v i (.elementSum (.getRow m_abs i))))
             (.maxElement v))
      :one (let [v ^AVector (Vector/createLength cols)
                 m_abs ^Matrix (m/abs m)]
             (dotimes [i cols]
               (.set v i (.elementSum (.getColumn m_abs i))))
             (.maxElement v))
      :two (norm m))))


(defn ^Double mat_complex_norm [^Matrix m & {:keys [type] :or {type :inf}}]
  (mat_norm (complexToReal m :type :mag) :type type))
