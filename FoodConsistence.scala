package AFSA

import scala.math._

/**
 * Created by jnu on 2015/12/22.
 */
class FoodConsistence {
  //适应度函数，求最大
  def function(x:Array[Double]):Double = {
    val m = x.length //变量x的个数
    var out: Double = 0
    out = -5*sin(x(0))*sin(x(1))*sin(x(2))*sin(x(3))*sin(x(4))-sin(5*x(0))*sin(5*x(1))*sin(5*x(2)) *sin(5*x(3))*sin(5*sin(4))+8
    1/out
  }
}
