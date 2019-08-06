#' Node class
#'
#' @field index numeric. Index of node in tree structure. Tree is indexed from top to bottom, left to right.
#' @field parent Node object. Points to parent node of given node.
#' @field leftChild Node object. Points to left child of given node.
#' @field rightChild Node object. Points to right child of given node.
#' @field weight numeric. Numeric weight of a given node. E.g. if tree used to
#'  bin points into set ranges, how many points are in the range spanned by that node.
#' @field depth numeric. Numeric depth of the node within tree structure.
#' @field range ANY. Numeric array c(min,max) indicating range represented by the node, where min is minimum of range and
#' max is maximum of range. 
#'
#' @return Node object.
#'
#' @examples
Node <- setRefClass(

 Class = 'Node',
  fields = list(
    index = 'numeric',
    parent = 'ANY',
    leftChild = 'ANY',
    rightChild = 'ANY',
    weight = 'numeric',
    depth = 'numeric',
    range = 'ANY'
  ))
 
 Node$methods(
   initialize = function(parent=NULL, depth=0, index=0, range=NULL){
     .self$index <- index
     .self$parent <- parent
     .self$depth <- depth
     .self$leftChild <- NULL
     .self$rightChild <- NULL
     .self$weight <- 0
     .self$range <- range
   }
 )

 Node$methods(
#' Increment weight of node by 1
#'
#' @return NULL
   incrementWeight = function(){
     .self$weight <- weight + 1
   }
 )
 
 Node$methods(
#' Add child node to .self
#'
#' @param isLeftChild Boolean, TRUE if child to add is left child, FALSE if child to add is right child.
   addChild=function(isLeftChild=TRUE){
     indx <- calculateIndex(isLeftChild)
     rng <- calculateRange(isLeftChild)
     dpth <- .self$depth + 1
     child = Node$new(parent=.self, depth=dpth, index=indx, range=rng)
     if(isLeftChild){
       .self$leftChild <- child
     }
     else{
       .self$rightChild <- child
     }
    }
 )
   
 Node$methods(
#' Calculate index of a node given its parent node.
#'
#' @param isLeftChild Boolean. TRUE if finding index of left child, FALSE if finding index of right child.
#'
#' @return Numeric index of child node.
   calculateIndex = function(isLeftChild){
     parentIndex <- .self$index
     if (isLeftChild){
       return(2*parentIndex)
     }
     else{
       return(2*parentIndex+1)
     }
   }
 )
  
 Node$methods(
#' Calculate range
#' 
#' Calculates range of child node given range of parent node. Returned range will be half of parent range;
#' if left child will be bottom half, if right child will be upper half of parent range.
#'
#' @param isLeftChild Boolean. TRUE if finding range of left child, FALSE if finding range of right child.
#'
#' @return Numeric of length 2; c(min,max) where min is minumum of child range and max is maximum of child range.
   calculateRange = function(isLeftChild){
     if (!is.null(.self$range)){
       parentMin <- .self$range[1]
       parentMax <- .self$range[2]
       midPoint <- parentMin + (parentMax-parentMin)/2
       
       if (isLeftChild){
         min <- parentMin
         max <- midPoint
       }
       else{
         min <- midPoint
         max <- parentMax
       }
       return(c(min,max))
     }
     else{
       warning(sprintf("Parent node with index %s has unspecified range. Defaulting to range NULL", .self$index))
       return(NULL)
     }
   }
 )


Tree <- setRefClass(
  Class = 'tree',
  fields = list(
    range = 'numeric',
    nLeaves = 'numeric',
    treeDepth = 'numeric',
    head = 'ANY'
  ))

Tree$methods(
  #' Binary Tree Creation
  #'
  #' Creates a weighted binary tree. Each node is a bucket over some data range, and its weight
  #' will be incremented by 1 for every data point that is within this range. This initialization
  #' only initializes the tree structure itself, given the total range \code{range} of the data that
  #' will be inserted and \code{nLeaves}, the number of leaf nodes to bucket the data in.
  #'
  #' E.g. given range \code{range = c(0,4)} and \code{nLeaves=2}, the tree will be of depth 2 and
  #' have 2 leaves, with one bucketing data between \eqn{[0,2)}
  #' and one bucketing data between \eqn{[2,4)}.
  #'
  #' If a parent node has a left child with range \eqn{[min1, max1)} and a right child with range
  #' \eqn{[min2, max2)}, then the parent node will have range \eqn{[min1, max2)}.
  #' As currently instantiated, nLeaves is required to be a power of 2 (i.e. the constructed tree is a
  #' perfect binary tree).
  #'
  #' @param range Range: An array of two numeric values that indicate the range of data that will be bucketed
  #'  in the tree structure. The first value should be the minimum (inclusive) and the second should be the
  #'  maximum (exclusive).
  #' @param nLeaves Number of leaves: A whole, positive number, indicating the number of leaves the tree should have.
  #'
  #' @return
  #' @export
  #'
  #' @examples
  initialize = function(range, nLeaves){
    checkNLeaves(nLeaves)
    checkRange(range)
    
    .self$range <- range
    .self$nLeaves <- nLeaves
    .self$treeDepth <- calculateDepth(.self$nLeaves)
    
    .self$head <- NULL
    
    buildTree()
  }
)

Tree$methods(
#' Checks if number of leaves is a power of two or not
#'
#' @param nLeaves Numeric value indicating the number of leaves in a tree object.
#'
#' @return Boolean value. TRUE if nLeaves is a power of two. Stops and returns error otherwise.
  checkNLeaves = function(nLeaves){
  if (log2(nLeaves)%%1 == 0){
    return(TRUE)
  }
  else{
    stop("Number of leaves specified is not a power of two")
  }
})

Tree$methods( 
#' Calculates depth of tree given number of leaves nLeaves.
#'
#' Assumes that tree is perfect; i.e. has leaves filled in to complete depth.
#' Assumes that a tree with only a single root node has depth 1. 
#' 
#' @param nLeaves Numeric value indicating the number of leaves in a tree object.
#'
#' @return Numeric value of depth of tree. 
  calculateDepth = function(nLeaves){
  return(log2(nLeaves)+1)
})

Tree$methods(
  #ToDo: write method that calculates size of subtree at node
  # should store this at nodes for ease
  # note this is trivial here since trees are perfect;
  # would otherwise need to be a recursive method.
  calculateSize = function(currNode=.self$head){
    print("blah")
  }
)

Tree$methods(
#' Insert a child node directly under a specified parent node. If
#' parent node has no children, will insert a left child. If there
#' already is a left child, a right child will be inserted where possible,
#' or an error is raised.
#'
#' @param parent Node object.
  insertChild = function(parent){
    if (is.null(parent$leftChild)){
      parent$addChild(isLeftChild=TRUE)
    }
    else if (is.null(parent$rightChild)) {
      parent$addChild(isLeftChild=FALSE)
    }
    else{
      stop("Parent node already has two children; cannot insert child.")
    }
  }
)

Tree$methods(
#' Insert two children beneath a specified parent node.
#'
#' @param parent Node object.
  insertChildren = function(parent){
    insertChild(parent)
    insertChild(parent)
  }
)

Tree$methods(
  #' Recursively build tree object up to .self$depth
  #'
  #' @param currNode Current node, used for recursive call. 
  buildTree = function(currNode=NULL){
    # If no currNode, this means we are beginning tree construction.
    if (is.null(currNode)){
      # Create head of tree and set currNode to head. 
      .self$head <- Node$new(depth=1, index=1, range=.self$range)
      currNode <- .self$head
    }
    # If tree depth greater than current node depth, insert children and recursively build the tree.
    if (.self$treeDepth > currNode$depth){
      insertChildren(currNode)
      buildTree(currNode=currNode$leftChild)
      buildTree(currNode=currNode$rightChild)
    }
    else{
      return(NULL)
    }
  }
)

Tree$methods(
  #ToDo: update this so it's just a formatting function so it can be easily called with callSuper
#' Simple print function to print tree structure for debugging.. Recurses through the nodes of the tree and prints 
#' information associated with each node.
#'
#' @param currNode Node object, defaults to the head of the tree. 
  printTree = function(currNode=.self$head){
    if (identical(currNode,.self$head)){
      print("")
      print("Index, Left Child,  Right Child, Range")
    }
    if (!is.null(currNode)){
      print(c(currNode$index, currNode$leftChild$index, currNode$rightChild$index, currNode$range))
      printTree(currNode=currNode$leftChild)
      printTree(currNode=currNode$rightChild)
    }
    else{
      return(NULL)
    }
  }
)

Tree$methods(
  #function to add values in toAdd to the nodes in the Tree object.
  #the nth item in itemsToAdd is added to the node with index n.
  add = function(itemsToAdd){
    #ToDo
  }
)

Tree$methods(
  addHelper = function(itemToAdd, targetNodeIndex, currNode){
    if (targetNodeIndex == currNode$index){
      currNode$weight <- currNode$weight + itemToAdd;
    }
    # else if (!is.null(currNode)){
    #   #determine direction to traverse tree
    #   if targetNodeIndex
    # }
    else{
      return(NULL)
    }
  }
)


publicTreeStatistic <- setRefClass(
  Class = 'publicTreeStatistic',
  contains = 'tree',
  fields = list(
    data = 'ANY'
  ))

publicTreeStatistic$methods(
  initialize = function(data, range, nLeaves){
    callSuper(range, nLeaves)
    .self$data = data
    binData()
  }
)

publicTreeStatistic$methods(
  binData = function(){
    for (x in data){
      binDataPoint(x, .self$head)
    }
  }
)

publicTreeStatistic$methods(
  binDataPoint = function(x, currNode){
    if (!is.null(currNode)){
      currNode$incrementWeight()
      #print(traverseLeft(x, currNode))
      if (traverseLeft(x, currNode)){
        binDataPoint(x, currNode$leftChild)
      }
      else{
        binDataPoint(x, currNode$rightChild)
      }
    }
    else{
      return(NULL)
    }
  }
)

publicTreeStatistic$methods(
  traverseLeft = function(x, currNode){
    if (!is.null(currNode$leftChild)){
      min <- currNode$leftChild$range[1]
      max <- currNode$leftChild$range[2]
      if (x >= min && x < max){
        return(TRUE)
      }
      else{
        return(FALSE)
      }
    }
    #if null, doesn't matter which direction traversal happens.
    return(TRUE) 
  }
)


publicTreeStatistic$methods(
  printTree = function(currNode=.self$head){
    if (identical(currNode,.self$head)){
      print("")
      print("Index, Left Child,  Right Child, Range, weight")
    }
    if (!is.null(currNode)){
      print(c(currNode$index, currNode$leftChild$index, currNode$rightChild$index, currNode$range, currNode$weight))
      printTree(currNode=currNode$leftChild)
      printTree(currNode=currNode$rightChild)
    }
    else{
      return(NULL)
    }
  }
)


# efficientTree <- setRefClass(
#   Class = 'efficientTree'
# )

# ToDo: add accuracy stuff, data imputation, percentiles, alpha
# ToDo: n.bins needs to be changed to nBins.
dpTreeStatistic <- setRefClass(
  Class = 'dpTreeStatistic',
  contains = 'mechanismLaplace',
  methods = list(
    initialize = function(variable, n.bins, rng, epsilon=NULL, impute.rng=NULL){
      .self$name <- 'differentially private binary tree'
      .self$mechanism <- 'mechanismLaplace'
      .self$variable <- variable
      .self$n.bins <- n.bins
      .self$rng <- rng
      .self$accuracy <- accuracy
    },
    release = function(data){
      x <- data[,variable]
    }
  )
)

#' Generalized Laplace mechanism
#' 
#' Laplace mechanism as applied to a generic object which has a specified
#' add() function that describes how the noise should be added. 
#'
#' @return
#' @export
#'
#' @examples
mechanismGenLaplace <- setRefClass(
  Class = 'GenLaplace',
  contains = 'mechanism'
)

mechanismGenLaplace$methods(
  evaluate = function(fun, x, sens, postFun, ...) {
    x <- censordata(x, .self$var.type, .self$rng, .self$bins)
    x <- fillMissing(x, .self$var.type, impute.rng=.self$rng, categories=.self$bins)
    fun.args <- getFuncArgs(fun, inputList=list(...), inputObject=.self)
    input.vals = c(list(x=x), fun.args)
    true.val <- do.call(fun, input.vals)
    scale <- sens / .self$epsilon
    release <- true.val$add(dpNoise(n=length(true.val), scale=scale, dist='laplace'))
    out <- list('release' = release)
    out <- postFun(out, ...)
    return(out)
  }
)







