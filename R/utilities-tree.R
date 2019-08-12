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
#' @name incrementWeight
#' 
#' @return NULL
   incrementWeight = function(){
     .self$weight <- weight + 1
   }
 )
 
 Node$methods(
#' Add child node to .self
#' @name addChild
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
#' @name calculateIndex
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
#' @name calculateRange
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
  #' Framework to creates a binary tree. Each node is a bucket over some data range, and its weight
  #' will be incremented by 1 for every data point that is within this range. This object
  #' only initializes the tree structure itself, given the total range \code{range} of the data that
  #' will be inserted and \code{nLeaves}, the number of leaf nodes to bucket the data in. The 
  #' actual bucketing of data occurs within the publicTreeStatistic object that inherits from this 
  #' basic Tree structure.
  #'
  #' E.g. given range \code{range = c(0,4)} and \code{nLeaves=2}, Tree$new(range,nLeaves) 
  #' will create a tree of depth 2 that will have 2 leaves, with one bucketing data between \eqn{[0,2)}
  #' and one bucketing data between \eqn{[2,4)}.
  #'
  #' If a parent node has a left child with range \eqn{[min1, max1)} and a right child with range
  #' \eqn{[min2, max2)}, then the parent node will have range \eqn{[min1, max2)}.
  #' As currently instantiated, nLeaves is required to be a power of 2 (i.e. the constructed tree is a
  #' perfect binary tree).
  #' @name initialize
  #' @param range Range: An array of two numeric values that indicate the range of data that will be bucketed
  #'  in the tree structure. The first value should be the minimum (inclusive) and the second should be the
  #'  maximum (exclusive).
  #' @param nLeaves Number of leaves: A whole, positive number, indicating the number of leaves the tree should have.
  #'
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
#' @name checkNLeaves
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
#' @name calculateDepth
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
#' @name insertChild
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
#' @name insertChildren
#' @param parent Node object.
  insertChildren = function(parent){
    insertChild(parent)
    insertChild(parent)
  }
)

Tree$methods(
  #' Recursively build tree object up to .self$depth
  #' @name buildTree
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
#' @name printTree
#' @param currNode Node object, defaults to the head of the tree. 
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

Tree$methods(
#' Add values to the weights of nodes in Tree object.
#' 
#' The nth item in itemsToAdd is added to the weight of the node with index n.
#' @name addItems
#' @param itemsToAdd Array of numeric values.
#'
#' @examples 
#' t <- Tree$new(c(0,4),2) # this will be a tree with 2 leaves and 3 total nodes including the root
#' t$addItems(c(1,3,7)) # this will add weight 1 to the root, weight 3 to the node with index 2, and weight 7 to the node with index 3.
#' t$printTree() 
  addItems = function(itemsToAdd){
    for(i in 1:length(itemsToAdd)){
      print(i)
      print(itemsToAdd[i])
      addItem(itemsToAdd[i],i)
    }
  }
)

Tree$methods(
#' Add numeric value itemtoAdd to node's weight. 
#' @name addItem
#' @param itemToAdd Numeric value
#' @param targetNodeIndex Integer, target index of the node whose weight is to be added to.
  addItem = function(itemToAdd, targetNodeIndex){
    currNode <- .self$head
    addHelper(itemToAdd, currNode, targetNodeIndex)
    }
)

Tree$methods(
  addHelper = function(itemToAdd, currNode, targetNodeIndex){
    #If target node is currentNode, add itemtoAdd to the currentNode
    if (targetNodeIndex == currNode$index){
      currNode$weight <- currNode$weight + itemToAdd;
    }
    #Otherwise, determine whether or not to traverse the tree to the left or right
    #and then recurse on that subtree.
    else{
      if (traverseLeft(currNode$index, targetNodeIndex)){
        currNode <- currNode$leftChild
      }
      else{
        currNode <- currNode$rightChild
      }
      addHelper(itemToAdd, currNode, targetNodeIndex)
    }
  }
)

Tree$methods(
#' Determine whether or not the targetNode is to the left of the
#' current node (currNode) in the tree, assuming the targetNode is beneath the
#' current node.
#' @name traverseLeft
#' @param currNodeIndex Positive integer, index of current node in tree.
#' @param targetNodeIndex Positive integer, index of target node in tree.
#'
#' @return Boolean. TRUE if targetNode is to the left of currNode, FALSE otherwise.
  traverseLeft = function(currNodeIndex, targetNodeIndex){
    isPositiveInteger <- (currNodeIndex%%1==0) && (targetNodeIndex%%1==0) && (currNodeIndex > 0) && (targetNodeIndex >0)
    isProperlyOrdered <- targetNodeIndex > currNodeIndex
    if (isPositiveInteger && isProperlyOrdered){
      #determine depth of targetNode and currNode (note that we start depth at 1 not 0, hence the addition of the 1)
      targetNodeDepth <- floor(log2(targetNodeIndex))+1
      currNodeDepth <- floor(log2(currNodeIndex))+1
      
      #difference in depth between targetNode and currNode
      depthDiff <- targetNodeDepth - currNodeDepth

      #the left-most index of the subtree that has currNode at its head and leaves at the depth of targetNode
      leftMostIndex <- currNodeIndex * 2**depthDiff

      #width of this subtree
      subtreeWidth <- 2**depthDiff
      
      #cut-off value for leaves of this subtree that are reachable by traversing from currNode left rather than right down the tree.
      cutOff <- leftMostIndex + subtreeWidth/2

      if(targetNodeIndex < cutOff){
        return(TRUE)
      }
      else{
        return(FALSE)
      }
    }
  else{
    if (!isPositiveInteger){
      stop("Inputs must be positive integers.")
    }
    if (!isProperlyOrdered){
      stop("currNodeIndex must be smaller than targetNodeIndex.")
    }
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
#' Public (i.e., not differentially private) tree statistic
#' @name initialize
#' @param data Array of numeric data values to bin.
#' @param range Array of numeric data of length 2; the range that the desired histogram will span.
#'  Array of form c(min, max), where the first value, min, is the minimum of the range, and the second value, max is the maximum of the range.
#' @param nLeaves Numeric, the number of leaves in the tree, i.e. the number of histogram bins at their most granular.
#'
#' By a 'tree statistic' here we simply mean a sort of glorified histogram: imagine you have a perfect tree, and each of the leaves represent a 
#' range [x,y]. You could then bucket data points by incrementing the weight of a leaf by 1 for every data point that lies within the range
#' represented by that bucket. If we just binned things based on the leaves of the tree, this is equivalent to a normal histogram, where the
#' weights of the leaves is the count of each histogram bucket.
#' 
#' However, we could also more generally assign ranges and weights to all nodes of the tree, by saying that nodes span the ranges of their combined
#' children and have weights equal to their combined children's weights. Formally, If [x1,y1] is the range of a node's left child and  [x2,y2] is 
#' the range of a node's right child, then say that node has range [x1,y2]. Similarly, if a node's left child has weight A and its right
#' child has weight B, then that node will have weight A+B. 
#' 
#' This 'tree' statistic gives us a way to calculate the CDF of a function with fewer queries. This is not so useful in the context of a public
#' statistic, where we know the bin counts with certainty. But in the differentially private setting, fewer queries corresponds to a reduction in
#' error in the calculation of the CDF.
#'
#'  The below example will create a tree with 2 leaves that range [0,4). The root node has range [0,4),
  #' its left child will have range [0,2), and its right child will have range [2,4). 
  #' The data, c(1,1,4) is then binned into each of these nodes. Since the head includes all points, 
  #' it will have weight 5. Its left child will have weight 4, since there are 4 points in the data that
  #' lie in the range [0,2), and its right child will have weight 1, since there is 1 point in the data
  #' that lies in the range [2,4). 
  #'  
  #' @examples
  #' data <- c(1,1,1,1,3)
  #' t <- publicTreeStatistic$new(data,c(0,4),2)
  initialize = function(data, range, nLeaves){
    callSuper(range, nLeaves) #calls the general Tree initialization function
    .self$data = data
    binData()
  }
)

publicTreeStatistic$methods(
#' Bin the data points in .self$data into the tree.
#' @name binData
  binData = function(){
    for (x in .self$data){
      binDataPoint(x, .self$head)
    }
  }
)

publicTreeStatistic$methods(
#' Helper function to bin a single data point into the tree.
#' Traverses the tree structure until the data point has been binned at all levels of the tree.
#' @name binDataPoint
#' @param x Numeric. Data point to be binned
#' @param currNode Node object. The current node to bin point into.
  binDataPoint = function(x, currNode){
    # If currNode is not NULL, continue.
    if (!is.null(currNode)){
      # Bin x into node
      currNode$incrementWeight()
      # Determine whether or not to traverse to left or right child of node
      if (traverseLeft(x, currNode)){
        binDataPoint(x, currNode$leftChild)
      }
      else{
        binDataPoint(x, currNode$rightChild)
      }
    }
    # If currNode is NULL, we've reached the bottom of the tree and can end.
    else{
      return(NULL)
    }
  }
)

publicTreeStatistic$methods(
#' Determine whether or not to traverse left down tree to bin data point x.
#' 
#' At node currNode, if x is in the range of the left child, should traverse left
#' down tree and increment weight of the left child. Otherwise, should traverse right.
#' @name traverseLeft
#' @param x Numeric. Data point to be binned.
#' @param currNode Node object. The current node to bin point into.
#' @return Boolean. TRUE if should traverse left, FALSE else.
#'
#' E.g. given currNode node with range [0,4), its left child will have range [0,2) and its
#' right child will have range [2,4). Then data point x = 1 should be binned into the left child so 
#' traverseLeft(1,currNode) will return TRUE. 
#' 
#' If instead x = 2, it should be binned into the right child, so traverseLeft(2,currNode) would 
#' return FALSE.
#' 
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