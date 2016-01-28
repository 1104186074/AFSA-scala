package AFSA

import scala.collection.mutable.ArrayBuffer
import scala.util.Random

/**
 * Created by jnu on 2015/12/22.
 */
class AFSA {
  def fun(fishnum:Int,trynumber:Int,maxgen:Int,visual:Double,delta:Double,step:Double,xminmax:Array[Double]): Array[Double] ={
   var matrix_X=init(fishnum,xminmax)
    val food=new Array[Double](fishnum)
    for (i<- 0 until fishnum){
      food(i)=(new FoodConsistence).function(matrix_X(i))
    }
    var bestX=matrix_X(0)
    var bestfood=food(0)
    for (i<-1 until fishnum){
      if (food(i)>bestfood){
        bestfood=food(i)
        bestX=matrix_X(i)
      }
    }
    var trace=new Array[Double](maxgen)
    var x1=new Array[Double](xminmax.length/2)
    var y1=0.0
    var x2=new Array[Double](xminmax.length/2)
    var y2=0.0
    var gen=0
    while (gen<maxgen){
      println(gen)
      for (i<-0 until fishnum){
        x1=swarm(matrix_X,i,visual,step,delta,trynumber,xminmax,food)
        y1=(new FoodConsistence).function(x1)
        x2=follow(matrix_X,i,visual,step,delta,trynumber,xminmax,food)
        y2=(new FoodConsistence).function(x2)
        if (y1>y2){
          matrix_X(i)=x1
          food(i)=y1
        }
        else{
          matrix_X(i)=x2
          food(i)=y2
        }
      }
      val index=(0 until fishnum).toArray
      val maxindex=maxIndex(index,food)
      if (bestfood<food(maxindex)&& (bestfood-food(maxindex))<10e-3){
        bestfood=food(maxindex)
        bestX=matrix_X(maxindex)
      }
      else{
        gen=maxgen
      }
      trace(gen)=bestfood
      gen=gen+1
    }
    bestX.foreach(println)
    trace
  }
//鱼群初始化
  def init(fishnum:Int,xminmax:Array[Double]):Array[Array[Double]]={
    val Dim=xminmax.length/2
    val matrix_X=Array.ofDim[Double](fishnum,Dim)
    for (i<-0 until fishnum){
      for (j<-0 until Dim){
        matrix_X(i)(j)=Random.nextDouble()*(xminmax(2*j+1)-xminmax(2*j))
      }
      boundaryRepair3(matrix_X(i),xminmax,Dim)
    }
    matrix_X
  }
  def boundaryRepair3(matrix_X:Array[Double],xminmax:Array[Double],Dim:Int):Array[Double]={
    for (i<-0 until Dim){
      if (matrix_X(i)<xminmax(2*i)){
        matrix_X(i)=xminmax(2*i)
      }
      if (matrix_X(i)>xminmax(2*i+1)){
        matrix_X(i)=xminmax(2*i+1)
      }
    }
    matrix_X
  }
  //觅食行为
  /**
   *
   * @param matrix_Xi 当前人工鱼的位置
   * @param i  当前人工鱼的编号
   * @param visual  视野
   * @param step   步长
   * @param trynumber   尝试次数
   * @param xminmax     解的边界大小
   * @param foodLast    上一次迭代各人工鱼的食物浓度
   * @return
   */
  def prey(matrix_Xi:Array[Double],i:Int,visual:Double,step:Double,trynumber:Int,
            xminmax:Array[Double],foodLast:Array[Double]):Array[Double]={
    var xNext=new Array[Double](matrix_Xi.length)
    val Foodi=foodLast(i)
    var matrix_Xj=new Array[Double](matrix_Xi.length)
    var Foodj=0.0
    var flag=true
    var k=0
    while (k<trynumber){
      matrix_Xj=matrix_Xi.map(_+(2*(new Random).nextDouble()-1)*visual)
      Foodj=(new FoodConsistence).function(matrix_Xj)
      val temp=sumAndsubtraction(matrix_Xj,matrix_Xi,false)
      if (Foodi<Foodj){
        xNext=sumAndsubtraction(matrix_Xi,temp.map(_*step*(new Random).nextDouble()*norm(temp)),true)
        xNext=boundaryRepair(xNext,xminmax,matrix_Xi.length)
        k=trynumber
        flag=false
      }
    }
    //随机行为
    if (flag){
      matrix_Xj=matrix_Xi.map(_+(2*(new Random).nextDouble()-1)*visual)
      xNext=matrix_Xj
      xNext=boundaryRepair(xNext,xminmax,matrix_Xi.length)
    }
    val FoodNext=(new FoodConsistence).function(xNext)
//    var result=ArrayBuffer[Double]()
//    result++=xNext
//    result+=FoodNext
    xNext
  }
//聚群行为
  /**
   *
   * @param matrix_X
   * @param i
   * @param visual
   * @param step
   * @param delta
   * @param trynumber
   * @param xminmax
   * @param foodLast
   * @return
   */
  def swarm(matrix_X:Array[Array[Double]],i:Int,visual:Double,step:Double,
             delta:Double,trynumber:Int,xminmax:Array[Double],foodLast:Array[Double]):Array[Double]={
  val D=dist(matrix_X(i),matrix_X)
  val index = D.map(_<visual)
  val dim=index.length
  var nf=0
  val indexNum=ArrayBuffer[Int]()
  for (j<-0 until dim) {
    if (index(j).equals(true) &&j!=i) {
      nf += 1
      indexNum+= j
    }
  }
  val Xc=new Array[Double](matrix_X(0).length)//中心位置
  var result=ArrayBuffer[Double]()
    var xNext=new Array[Double](matrix_X(0).length)
  if (nf > 1){
    for (k<-0 until matrix_X(0).length){
      Xc(k)=mean(matrix_X,indexNum,k)
    }
    val Yc=(new FoodConsistence).function(Xc)
    val Yi=foodLast(i)

    if (Yc/(nf-1)>delta*Yi){
      val temp=sumAndsubtraction(Xc,matrix_X(i),false)
      xNext=sumAndsubtraction(matrix_X(i),temp.map(_*step*(new Random).nextDouble()*norm(temp)),true)
      xNext=boundaryRepair(xNext,xminmax,matrix_X(i).length)
//      val yNext=(new FoodConsistence).function(xNext)
//      result++=xNext
//      result+=yNext
    }
    else{
      xNext=prey(matrix_X(i),i,visual,step,trynumber,xminmax,foodLast)
    }
  }
  else{
    xNext=prey(matrix_X(i),i,visual,step,trynumber,xminmax,foodLast)
  }
  //result.toArray
    xNext
}
//追尾行为
  /**
   *
   * @param matrix_X
   * @param i
   * @param visual
   * @param step
   * @param delta
   * @param trynumber
   * @param xminmax
   */
  def follow(matrix_X:Array[Array[Double]],i:Int,visual:Double,step:Double,delta:Double,
              trynumber:Int,xminmax:Array[Double],food:Array[Double]):Array[Double]={
    val D=dist(matrix_X(i),matrix_X)
    val index = D.map(_<visual)
    val dim=index.length
    var nf=0
    val indexNum=ArrayBuffer[Int]()
    for (j<-0 until dim) {
      if (index(j).equals(true)&&j!=i)
        nf += 1
      indexNum+=j
    }
    var result=ArrayBuffer[Double]()
    var xNext=new Array[Double](dim)
    if (nf > 1){
      val Xmax=matrix_X(maxIndex(indexNum.toArray,food))
      val Ymax=food(maxIndex(indexNum.toArray,food))
      val Yi=food(i)
      if (Ymax/(nf-1)>delta*Yi){
        val temp=sumAndsubtraction(Xmax,matrix_X(i),false)
        xNext=sumAndsubtraction(matrix_X(i),temp.map(_*step*(new Random).nextDouble()*norm(temp)),true)
        xNext=boundaryRepair(xNext,xminmax,matrix_X(i).length)
//        val yNext=(new FoodConsistence).function(xNext)
//        result++=xNext
//        result+=yNext
      }
      else{
        xNext=prey(matrix_X(i),i,visual,step,trynumber,xminmax,food)
      }
    }
    else{
      xNext=prey(matrix_X(i),i,visual,step,trynumber,xminmax,food)
    }
    xNext
    //result.toArray
  }


  //矩阵相加减法
  /**
   *
   * @param array1 矩阵1
   * @param array2 矩阵2
   * @param b true为加法，false为减法
   * @return
   */
  def sumAndsubtraction(array1: Array[Double],array2: Array[Double],b: Boolean):Array[Double]={
    val result=new Array[Double](array1.length)
    if (array1.length==array2.length) {
      if (b.equals(true)) {
        for (i <- 0 until array1.length) {
          result(i) = array1(i) + array2(i)
        }
      }
      else {
        for (i <- 0 until array1.length) {
          result(i) = array1(i) - array2(i)
        }
      }
    }
    result
  }

  //二范数
  def norm(Array:Array[Double]):Double={
    var result=0.0
    for (i<-0 until Array.length){
      result+=math.pow(Array(i),2)
    }
    math.sqrt(result)
  }
//边界修正
  def boundaryRepair(matrix_X:Array[Double],xminmax:Array[Double],Dim:Int):Array[Double]={
    for (i<-0 until Dim){
      if (matrix_X(i)<xminmax(2*i)){
        matrix_X(i)=xminmax(2*i)
      }
      if (matrix_X(i)>xminmax(2*i+1)){
        matrix_X(i)=xminmax(2*i+1)
      }
    }
    matrix_X
  }

  //计算一条鱼与所有鱼的距离
  /**
   *
   * @param matrix_Xi  第i条鱼的当前位置
   * @param matrix_X   所有鱼的当前位置
   * @return   距离数组
   */
  def dist(matrix_Xi:Array[Double],matrix_X:Array[Array[Double]]):Array[Double]={
    val popsize=matrix_X.length
    val distance=new Array[Double](popsize)
    for (i<-0 until popsize){
      distance(i)=norm(sumAndsubtraction(matrix_Xi,matrix_X(i),false))
    }
    distance
  }

  //求平均值
  def mean(matrix_X:Array[Array[Double]],indexNum:ArrayBuffer[Int],k:Int):Double={
    var result=0.0
    for (i<- 0 until indexNum.length){
      result+=matrix_X(indexNum(i))(k)
    }
    result/(indexNum.length)
  }
  //求最大食物浓度对应的index
  def maxIndex(indexNum:Array[Int],food:Array[Double]):Int={
    var maxindex=indexNum(0)
    var foodmax=food(maxindex)
    for (i<- 1 until indexNum.length){
      if (foodmax<food(indexNum(i))){
        maxindex=indexNum(i)
        foodmax=food(maxindex)
      }
    }
    maxindex
  }

}

