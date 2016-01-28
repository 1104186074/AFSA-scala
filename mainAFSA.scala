package AFSA

/**
 * Created by jnu on 2015/12/23.
 */
import scala.math._
object mainAFSA {
  def main(args: Array[String]) {
    val fishnum = 50
    val maxgen = 20
    val trynumber = 5
    val visual = 2.5
    val delta = 0.618
    val step = 0.3
    val xminmax = Array[Double](0, 0.9 * Pi, 0, 0.9 * Pi, 0, 0.9 * Pi, 0, 0.9 * Pi, 0, 0.9 * Pi)
    val result = (new AFSA).fun(fishnum, trynumber, maxgen, visual, delta, step, xminmax)
    result.map(1 / _).foreach(println)
  }
}
