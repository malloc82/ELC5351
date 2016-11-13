(ns elc5351.fft
  (:use clojure.core)
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.operators :as op]
            [clojure.core.matrix.linear    :refer [norm]])
  (:import [mikera.matrixx Matrix Matrixx]
           [mikera.vectorz AVector Vector Vectorz]
           [org.jtransforms.fft DoubleFFT_1D DoubleFFT_2D]))

(set! *warn-on-reflection* true)

(m/set-current-implementation :vectorz)

(defn ^Matrix createComplexRandom
  ([^Integer rows]
   (createComplexRandom rows rows))
  ([^Integer rows ^Integer cols]
   (Matrix/createRandom  rows (bit-shift-left cols 1))))

(defn ^Matrix realToComplex [^Matrix m]
  (let [[rows cols] (m/shape m)
        m_complex  ^Matrix  (Matrix/create rows (bit-shift-left cols 1))
        m_arr      ^doubles (.asDoubleArray m)
        m_comp_arr ^doubles (.asDoubleArray m_complex)]
    (dotimes [i (alength m_arr)]
      (aset m_comp_arr (bit-shift-left i 1) (aget m_arr i)))
    m_complex))

(defn ^Matrix complexToReal [^Matrix m & {:keys [type]
                                          :or {type :mag}}]
  (let [[rows cols] (m/shape m)
        cols_real   ^Integer (bit-shift-right cols 1)
        real_mat    ^Matrix  (Matrix/create rows cols_real)
        real_arr    ^doubles (.asDoubleArray real_mat)
        complex_arr ^doubles (.asDoubleArray m)]
    (case type
      :mag   (dotimes [i (* rows cols_real)]
               (let [idx (bit-shift-left i 1)]
                 (aset real_arr i ^Double (norm [(aget complex_arr idx)
                                                 (aget complex_arr (unchecked-inc idx))]))))
      :phase (dotimes [i (* rows cols_real)]
               (let [idx (bit-shift-left i 1)]
                 (aset real_arr i ^Double (Math/atan (/ (aget complex_arr (unchecked-inc idx))
                                                        (aget complex_arr idx)))))))
    real_mat))

(defn fft2 [^Matrix m & {:keys [complex copy]
                         :or {copy true}}]
  (let [[rows cols]         (m/shape m)
        complex_cols        (if complex (bit-shift-right cols 1) cols)
        complex_mat ^Matrix (if complex
                              (if copy (.clone m) m)
                              (realToComplex m))]
    (m/pm complex_mat)
    (.complexForward (DoubleFFT_2D. rows complex_cols) (.asDoubleArray complex_mat))
    complex_mat))


(defn ifft2 [^Matrix m & {:keys [copy scaling]
                          :or {copy true
                               scaling true}}]
  (let [[rows cols] (m/shape m)
        out (if copy (.clone m) m)]
    (.complexInverse ^DoubleFFT_2D (DoubleFFT_2D. rows (bit-shift-right cols 1))
                     (.asDoubleArray out) true)
    out))
