Node <- setRefClass(

 Class = 'Node',
  fields = list(
    index = 'numeric',
    parent = 'ANY',
    leftChild = 'ANY',
    rightChild = 'ANY',
    weight = 'numeric',
    range = 'ANY'
  ),

 methods = list(
   initialize = function(parent=NULL, index=0, range=NULL){
     .self$index <- index
     .self$parent <- parent
     .self$leftChild <- NULL
     .self$rightChild <- NULL
     .self$weight <- 0
     .self$range <- range
   },
   
   incrementWeight = function(){
     .self$weight <- weight + 1
   },
   addLeftChild = function(){
     indx <- calculateIndex(isLeftChild=TRUE)
     rng <- calculateRange(isLeftChild=TRUE)
     child <- Node$new(parent=.self, index=indx, range=rng)
     .self$leftChild = child
   },
   addRightChild = function(){
     indx <- calculateIndex(isLeftChild=FALSE)
     rng <- calculateRange(isLeftChild=FALSE)
     child = Node$new(parent=.self, index=indx, range=rng)
     .self$rightChild <- child
   },
   calculateIndex = function(isLeftChild){
     parentIndex <- .self$index
     if (isLeftChild){
       return(2*parentIndex)
     }
     else{
       return(2*parentIndex+1)
     }
   },
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
)

Tree <- setRefClass(
  Class = 'tree',
  fields = list(
    range = 'numeric',
    nLeaves = 'numeric',
    depth = 'numeric',
    head = 'ANY'
  ),
  methods = list(
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
    
    insertChild = function(parent){
      if (is.null(parent$leftChild)){
        parent$addLeftChild()
      }
      else if (is.null(parent$rightChild)) {
        parent$addRightChild()
      }
      else{
        stop("Parent node already has two children; cannot insert child.")
      }
    },
    
    insertChildren = function(parent){
      insertChild(parent)
      insertChild(parent)
    },
    #ToDo: update so that depth isn't passed since it's unnecessary. can just make sure currNode not null.
    buildTree = function(currDepth=1, currNode=NULL){
      if (currDepth == 1){
        .self$head <- Node$new(index=1, range=.self$range)
        currNode <- .self$head
      }
      if (.self$depth > currDepth){
        insertChildren(currNode)
        buildTree(currDepth+1, currNode=currNode$leftChild)
        buildTree(currDepth+1, currNode=currNode$rightChild)
      }
      else{
        return(NULL)
      }
    },
    
    #update this so it's just a formatting function so it can be easily called with callSuper
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
  )

publicTreeStatistic <- setRefClass(
  Class = 'publicTreeStatistic',
  contains = 'tree',
  fields = list(
    data = 'ANY'
  ),
  methods = list(
    initialize = function(data, range, nLeaves){
      callSuper(range, nLeaves)
      .self$data = data
      binData()
    },
    binData = function(){
      for (x in data){
        binDataPoint(x, .self$head)
      }
    },
    binDataPoint = function(x, currNode){
      if (!is.null(currNode)){
        currNode$incrementWeight()
        print(traverseLeft(x, currNode))
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
    },
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
    },
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
)





