Node <- setRefClass(
  
 Class = 'Node',
  fields = list(
    parent = 'ANY',
    leftChild = 'ANY',
    rightChild = 'ANY',
    weight = 'numeric',
    rng = 'ANY'
  ),
 
 methods = list(
   initialize = function(parent=NA){
     .self$parent <- parent
     .self$leftChild <- NA
     .self$rightChild <- NA
     .self$weight <- 0
     .self$rng <- NULL
   },
   setWeight = function(weight){
     .self$weight = weight
   },
   addLeftChild = function(){
     child = Node$new(parent=.self)
     .self$leftChild = child
   },
   addRightChild = function(){
     child = Node$new(parent=.self)
     .self$rightChild = child
   },
 )
)

Tree <- setRefClass(
  Class = 'tree',
  fields = list(
    rng = 'numeric',
    nLeaves = 'numeric',
    depth = 'numeric',
    head = 'ANY'
  ),
  methods = list(
    #' Binary Tree Creation
    #' 
    #' Creates a weighted binary tree. Each node is a bucket over some data range, and its weight
    #' will be incremented by 1 for every data point that is within this range. This initialization
    #' only initializes the tree structure itself, given the total range \code{rng} of the data that
    #' will be inserted and \code{nLeaves}, the number of leaf nodes to bucket the data in.
    #'
    #' E.g. given range \code{rng = c(0,4)} and \code{nLeaves=2}, the tree will be of depth 2 and
    #' have 2 leaves, with one bucketing data between \eqn{[0,2)} 
    #' and one bucketing data between \eqn{[2,4)}. 
    #' 
    #' If a parent node has a left child with range \eqn{[min1, max1)} and a right child with range
    #' \eqn{[min2, max2)}, then the parent node will have range \eqn{[min1, max2)}.
    #' As currently instantiated, nLeaves is required to be a power of 2 (i.e. the constructed tree is a
    #' perfect binary tree). 
    #' 
    #' @param rng Range: An array of two numeric values that indicate the range of data that will be bucketed
    #'  in the tree structure. The first value should be the minimum (inclusive) and the second should be the
    #'  maximum (exclusive).
    #' @param nLeaves Number of leaves: A whole, positive number, indicating the number of leaves the tree should have.
    #' 
    #' @return
    #' @export
    #'
    #' @examples
    initialize = function(rng, nLeaves){
      checkNLeaves(nLeaves)
      checkRange(rng)
      
      .self$rng <- rng
      .self$nLeaves <- nLeaves
      .self$depth <- calculateDepth(.self$nLeaves)
      .self$head <- NULL
      
      buildTree()
    },
    
    
    checkNLeaves = function(nLeaves){
      if (log2(nLeaves)%%1 == 0){
        return(TRUE)
      }
      else{
        stop("Number of leaves specified is not a power of two")
      }
    },
    
    calculateDepth = function(nLeaves){
      return(log2(nLeaves)+1)
    },
    
    buildTree = function(currDepth=1, currNode=NULL){
      if (currDepth == 1){
        .self$head <- Node$new()
        currNode <- .self$head
      }
      if (.self$depth > currDepth){
        currNode$addLeftChild()
        currNode$addRightChild()
        buildTree(currDepth+1, currNode=.self$head$leftChild)
        buildTree(currDepth+1, currNode=.self$head$rightChild)
      }
      else{
        return(NULL)
      }
    },
    
    printTree = function(currDepth=1, currNode=.self$head){
      if(currDepth=1){
        print()
        print($self)
      }
    }
  )
  )



