# John Nico T. De Castro
# 2022-12523
# B1L     Project
# Generic Solvers

# ----- Main Code for PolynomialRegression() -----


# PolynomialRegression()
## requirements: integer specifying order of polynomial and list of two vectors for datapoints
## results:
##      checks initial conditions for order and the length of dependent/independent variables
##      proceeds calling helper functions for the vandermonde matrix, gauss jordan method, and polynomial string and function
##      returns a list variable containing the augcoeffmatrix (vandermondematrix), coefficients, and the polynomial string and function
PolynomialRegression <- function(order, dataPoints, x) {
  n = length(dataPoints[[1]])
  if(order < 1 || order >= n) return(NA) # check if order is invalid (must be 1 < order < n)
  else if(n != length(dataPoints[[2]])) return(NA)  # check length of dependent/independent variables are different
  
  vandermondematrix = VandermondeMatrix(order, dataPoints)  # given order and datapoints generate the vandermondematrix
  
  gaussjordanresult = GaussJordanMethod(vandermondematrix)  # given the vandermondematrix, perform gauss jordan method
  
  polynomial_string = convertPolynomialString(gaussjordanresult$solution) # given the solution to gauss jordan, convert to polynomial string
  
  polynomial_function = eval(parse(text=polynomial_string)) # convert polynomial string to function
  
  return(list(coefficients = gaussjordanresult$solution, function1 = polynomial_function, estimate = polynomial_function(x), fxn_string = polynomial_string))
}


# VandermondeMatrix()
## requirements: integer specifying order of polynomial and list of two vectors for datapoints
## results:
##      set m and set empty vandermonde matrix
##      builds the vandermonde matrix horizontally using a helper function, summationX()
##      generates the variables using another helper function, generateVariables()
##      returns a list containing the variables and augcoeffmatrix
VandermondeMatrix <- function(order, dataPoints) {
  m = order+1
  vandermondematrix = matrix(NA, nrow=m, ncol=m+1, byrow=TRUE)
  
  # build vandermonde matrix horizontally
  for(row in 0:(order)) {   # compute for the summation of x^n and assigns it to the row of the matrix
    vandermondematrix[row+1,] = summationX(row, order+row, dataPoints)
  }
  
  # generate variables
  vars = generateVariables(order+1)
  
  return(list(variables = vars, augcoeffmatrix = vandermondematrix))
}


# summationX()
## requirements: current row, end exponent for the current row, and the datapoints
## results:
##      computes the lhs of the current row (summation of ind var raised by expo)
##      computes the rhs of the current row (summation of ind var raised by expo multiplied by y)
##      returns the computed row
summationX <- function(currRow, endExponent, dataPoints) {
  row = c()
  for(expo in currRow:endExponent) {  # iterates from current row to end exponent, and computes for the sum of the ind var raised by expo
    row = c(row, sum(dataPoints[[1]]**expo))
  }
  row = c(row, summationXY(currRow, dataPoints))  # computation for rhs, summation of ind var raised to the row no. multiplied by dep var
  return(row)
}


# summationXY()
## requirements: current row and the datapoints
## results:
##      returns the computed rhs of the current row
summationXY <- function(row, dataPoints) {
  x = dataPoints[[1]]
  y = dataPoints[[2]]
  sumXY = sum((x**row)*y)
  return(sumXY)
}


# generateVariables()
## requirements: order or degree of the polynomial
## results:
##      returns the generated vector of variables of the vandermonde matrix (requisite for gaussjordanmethod)
generateVariables <- function(order) {
  vars = c()
  for(var in 1:order) {
    vars = c(vars, paste("x", var, sep=""))
  }
  return(vars)
}


# convertPolynomialString()
## requirements: solutions vector
## results:
##      adds first solution to polynomial vector
##      calls the helper function, generatePolynomialVariables to generate the variables
##      iteratively adds each solution in solutions vector and the associated variable to the poly string vector (also converts it to a string vector) 
##      returns the polynomial in string form
convertPolynomialString <- function(solution) {
  polynomial_string = c(paste(solution[1], "*x^0", sep="", collapse=""))
  noOfVars = length(solution)
  vars = generatePolynomialVariables(noOfVars)
  
  for(currVar in 2:noOfVars) {    # Ex: 2nd term in solution and x ^ 1, 3rd term in solution and x ^ 2, ...
    polynomial_string = c(polynomial_string, paste(solution[currVar], vars[currVar-1], sep="*", collapse=""))
    
  }
  
  polynomial_string = paste(polynomial_string, collapse = " + ")    # separate each term with +
  
  return(paste("function(x)", polynomial_string, collapse = ""))    # concatenates the function(x) string
}


# generateVariables()
## requirements: noOfVars to be generated
## results:
##      returns the generated string vector of variables of the polynomial (in form x ^ expo)
generatePolynomialVariables <- function(noOfVars) {
  vars = c()
  for(exp in 1:noOfVars) {
    vars = c(vars, paste("x^", exp, sep="", collapse=""))
  }
  return(vars)
}


# ----- Main Code for GaussJordanMethod() ----- #

GaussJordanMethod <- function(augCoeffList) {
  if(anyNA(augCoeffList)) return(NA)    # if augCoeffList returned NA (not square matrix or # of frequency of unique vars are not equal)
  n = length(augCoeffList$variables)
  coeffMatrix = augCoeffList$augcoeffmatrix
  for(i in 1:n) {
    if(i != n) {
      # find index of pivot row (row with max(abs(a[i:n, i]))) using findPivotRowIndex helper function
      pivotRowIndex = findPivotRowIndex(coeffMatrix, i, n)
      
      # if no unique solution
      if(coeffMatrix[pivotRowIndex, i] == 0) return(NA)
      
      # partial pivoting, using swap helper function if current pivot row is not the correct pivot row
      if(i != pivotRowIndex) coeffMatrix = swap(coeffMatrix, i, pivotRowIndex)
    }
    coeffMatrix[i, ] = coeffMatrix[i, ] / coeffMatrix[i, i]     # normalize pivot row (pivot row/pivot element)
    for(j in 1:n) {   
      if(i == j) next
      normalizedRow = coeffMatrix[j, i] * coeffMatrix[i, ]  # get normalized row (value to elim * pivot row)
      coeffMatrix[j, ] = coeffMatrix[j, ] - normalizedRow  # original row - normalized row
    }
    # check main diagonal after normalizing column
    if(!coeffMatrix[i, i] %in% 1) return(NA)      # if main diagonal is not == 1, use %in% to handle NaN values
  }
  
  # build solution vector
  solution = c() # create empty solutions vector
  for(i in 1:n) {
    solution = c(solution, coeffMatrix[i, n+1])
  }
  
  return(list(variables = augCoeffList$variables, augcoeffmatrix = coeffMatrix, solution = solution))
}


# findPivotRowIndex()
## requirements: the coeffMatrix, current i, and the number of rows
## results:
##      gets absolute value of the ith column
##      orders the coefficients in non-descending order in a list of indices (rowIndex)
##      returns index of the pivot row, (1st in list + the offset)
findPivotRowIndex <- function(coeffMatrix, i, n) {
  coeffMatrix = abs(coeffMatrix[i:n, i])        # get absolute value of ith column of coeffMatrix
  rowIndex = order(coeffMatrix, decreasing = TRUE)           # get indices of each row in non-descending order
  return(rowIndex[1]+(i-1))           # get the index of row with max value (1st since non-descending order), offset by i-1 since 2nd iteration starts at 2, 3rd starts at 3
}


# swap()
## requirements: the coeffMatrix, the current row index, and correct pivot row index
## results:
##      performs simple swapping using a temp vector
##      returns the coeffMatrix with swapped rows
swap <- function(coeffMatrix, currentRowIndex, pivotRowIndex) {
  tempRow = coeffMatrix[currentRowIndex, ]
  coeffMatrix[currentRowIndex, ] = coeffMatrix[pivotRowIndex, ]
  coeffMatrix[pivotRowIndex, ] = tempRow
  return(coeffMatrix)
}


# ----- Main Code for AugCoeffMatrix() ----- #
# AugCoeffMatrix()
## requirements: list containing the equations following the rules
## results:
##      checks first if the # of unknown variables are equal
##      builds the matrix by row using extractCoeff() on each equation
##      transpose RHS
##      returns the variables involved and augmented coefficient matrix in a single list variable 
AugCoeffMatrix <- function(system) {
  if(isFreqNotEqual(system)) return(NA);    # checks if frequency of unique variables are not equal using isFreqNotEqual(), returns NA if TRUE
  vars = extractVar(system[[1]])      # extracts variables using extractVar() on first equation
  noOfFuncs = length(system)
  noOfVars = length(vars)
  if(noOfFuncs != noOfVars) return(NA) # if results to a non-square matrix, return na
  buildmatrix = matrix(NA, nrow=noOfFuncs, ncol=noOfVars+1, byrow=TRUE, dimnames=list(1:noOfFuncs,c(vars,"RHS")))
  
  for(row in 1:noOfFuncs) {     # builds the matrix by row, calling extractCoeff() on each equation
    buildmatrix[row,] = extractCoeff(system[[row]], noOfVars)
  }
  
  # transpose RHS 
  for(row in 1:nrow(buildmatrix)) {   # iterates over number of rows
    buildmatrix[row,(noOfVars+1)] = -buildmatrix[row,(noOfVars+1)]    # assign negated values in last column
  }
  
  return(list(variables = vars, augcoeffmatrix = buildmatrix))
}


# isFreqNotEqual()
## requirements: list containing the linear eqs
## results:
##      calls helper function getFrequency()
##      compares the frequency of the first variable to every other variable
##      returns TRUE if a different frequency is found, otherwise returns FALSE
isFreqNotEqual <- function(system) {
  matrixFrequency = getFrequency(system)
  for(i in 1:nrow(matrixFrequency)) {
    if(matrixFrequency[1,2] != matrixFrequency[i,2]) return(TRUE)
  }
  return(FALSE)
}


# getFrequency()
## requirements: list containing the equations
## results:
##      extracts variables from each equation using extractVar() and assigns it to a vector, varList
##      checks frequency of each unique variable
##      returns matrixFrequency (frequencies of variables in matrix form)
getFrequency <- function(system) { # modified getFrequency algo from prev. exer
  uniqueVars = c()
  varFrequencies = numeric()
  for(i in 1:length(system)) {
    varList = extractVar(system[[i]])
    for(var in varList) {
      if(var %in% uniqueVars) {
        varIndex = which(uniqueVars==var)
        varFrequencies[varIndex] = varFrequencies[varIndex] + 1
      } else {
        uniqueVars = c(uniqueVars, var)
        varFrequencies = c(varFrequencies, 1)
      }
    }
  }
  
  matrixFrequency <- matrix(
    c(uniqueVars, varFrequencies),
    nrow = length(uniqueVars),
    ncol = 2,
    byrow = FALSE,
    dimname = list(1:length(uniqueVars),c("Unique Variable", "Frequency"))
  )
  return(matrixFrequency)
}


# extractVar()
## requirements: a single function
## results:
##      extracts the variables from the function parameters
##      returns these variables in a character vector, varList
extractVar <- function(func) {
  varList = deparse(func)[1]      # deparses func (expression into character strings) and returns the 1st element: function (x1, x2, ..., xn)
  varList = substring(varList, 10)      # removes function, returns only the opening parenthesis until end of string
  varList = gsub("[() ]", "", varList)      # removes spaces and parentheses
  varList = strsplit(varList, ",")[[1]]     # splits using comma as delimiter, [[1]] required since it splits only one string
  return(varList)
}


# extractCoeff()
## requirements: a single function and number of variables
## results:
##      calls helper function extractedSortedItems()
##      extracts the coefficients from termList
##      returns coeffList (numeric vector of coefficients)
extractCoeff <- function(func, noOfVars) {
  termList = extractSortedTerms(func, noOfVars) # helper func that returns termList (vector of terms in the expression)
  coeffList = c()
  for(term in termList) {     # iterates over terms, extracts coefficient, and adds it to coeffList
    coeffList = c(coeffList, gsub("\\*x[0-99]", "", term))
  }
  return(as.numeric(coeffList))
}


# extractedSortedTerms()
## requirements: a single function and number of variables
## results:
##      deparses an input function, removes spaces and splits by the operator +
##      sorts the character using a selection sort algorithm
##      returns Termlist (vector of terms in the expression)
extractSortedTerms <- function(func, noOfVars) {
  termList = deparse(func, width.cutoff = 200)[2]   
  termList = gsub(" ", "", termList)
  termList = strsplit(termList, "\\+")[[1]]
  for(i in 1:(noOfVars-1)) {      # selection sort algo
    minIndex = i
    for(j in (1+i):(noOfVars)) {
      currentTerm = as.integer(substring(termList[minIndex],nchar(termList[minIndex]))[[1]])
      nextTerm = as.integer(substring(termList[j],nchar(termList[j]))[[1]])
      if(currentTerm > nextTerm) {
        minIndex = j
      }
    }
    if(minIndex != i) {     # swap positions in termList
      temp = termList[minIndex]
      termList[minIndex] = termList[i]
      termList[i] = temp;
      
    }
  }
  return(termList)
}


# ----- Main Code for QSI() ----- #


QSI <- function(mat, x) {
  mat = mat[order(mat[,1]),]
  n = nrow(mat)-1     # number of intervals
  numOfEqs = 3*n      # number of total num. of eqs
  
  eqsMatrix = matrix(0, nrow=numOfEqs-1, ncol=numOfEqs)
  
  # condition 1
  currEq = 1
  for(i in 2:n) {
    b = mat[i, 1]
    a = b**2
    rhs = mat[i, 2]
    
    # a
    if(i-1 != 1) {
      eqsMatrix[currEq,i-2] = a
    }
    eqsMatrix[currEq+1,i-1] = a
    
    # b
    eqsMatrix[currEq, i+(n-2)] = b
    eqsMatrix[currEq+1, i+(n-1)] = b
    
    # c
    eqsMatrix[currEq, i+(n*2-2)] = 1
    eqsMatrix[currEq+1, i+(n*2-1)] = 1
    
    # rhs
    eqsMatrix[currEq, numOfEqs] = rhs
    eqsMatrix[currEq+1, numOfEqs] = rhs
    
    currEq = currEq + 2
  }
  
  # condition 2
  # a
  eqsMatrix[currEq+1, n-1] = mat[n+1, 1]**2
  
  # b
  eqsMatrix[currEq, n] = mat[1, 1]
  eqsMatrix[currEq+1, (n*2)-1] = mat[n+1, 1]
  
  # c
  eqsMatrix[currEq, n*2] = 1
  eqsMatrix[currEq+1, (n*3)-1] = 1
  
  # rhs
  eqsMatrix[currEq, numOfEqs] = mat[1, 2]
  eqsMatrix[currEq+1, numOfEqs] = mat[n+1, 2]
  currEq = currEq + 2
  
  # condition 3
  for(i in 2:n) {
    # a
    a = mat[i,1]*2
    if(i-1 != 1) {
      eqsMatrix[currEq, i-2] = a
    }
    eqsMatrix[currEq, i-1] = -a
    
    # b
    eqsMatrix[currEq, i+(n-2)] = 1
    eqsMatrix[currEq, i+(n-1)] = -1
    
    currEq = currEq + 1
  }

  
  # gaussjordanmethod 
  gaussJordanResult = ModifiedGaussJordanMethod(eqsMatrix)
  
  # build intervals
  intervals = buildIntervals(gaussJordanResult$solution, n)
  
  # choose the appropriate interval
  if(x > max(mat[,1]) || x < min(mat[,1])) {
    x_interval = "x value is not within the given intervals!"
    fxn = "x value is not within the given intervals!"
    fxn_string = "x value is not within the given intervals!"
    estimate = "x value is not within the given intervals!"
  } else {
    for(i in 1:n) {
      if(mat[i, 1] <= x && x <= mat[i+1, 1]) {
        x_interval = paste(mat[i,1], "<=", "x", "<=", mat[i+1,1])
        fxn_string = intervals[i]
        break
      }
    }
    
    # build usable function
    fxn = eval(parse(text=fxn_string))
    
    #estimate
    estimate = fxn(x)
  }
  
  
  return (list(intervals=intervals,fxn=fxn,estimate=estimate,fxn_string=fxn_string,x_interval=x_interval))
}


buildIntervals <- function(solution, n) {
  # 0 3 6
  # 1 4 7
  # 2 5 8
  intervals = c()
  curr = 0
  coeffVector = c()
  for(i in 1:n) {
    if(curr != 0) {
      coeffVector = c(coeffVector, paste(solution[curr], "*(x^2)", sep=""))
    }
    coeffVector = c(coeffVector, paste(solution[curr+n], "*(x^1)", sep=""))
    coeffVector = c(coeffVector, paste(solution[curr+n*2], "", sep=""))
    curr = curr+1
  }
  intervals = c(intervals, paste("function(x)", paste(coeffVector[1], coeffVector[2], sep=" + ", collapse=""), collapse=""))
  counter = 3
  for(i in 2:n) {
    intervals = c(intervals, paste("function(x)", paste(coeffVector[counter], coeffVector[counter+1], coeffVector[counter+2], sep=" + ", collapse=""), collapse=""))
    counter = counter + 3
  }
  return(intervals)
}


# ----- Main Code for GaussJordanMethod() ----- #

ModifiedGaussJordanMethod <- function(augCoeffMatrix) {
  if(anyNA(augCoeffMatrix)) return(NA)    # if augCoeffList returned NA (not square matrix or # of frequency of unique vars are not equal)
  n = ncol(augCoeffMatrix)-1
  coeffMatrix = augCoeffMatrix
  for(i in 1:n) {
    if(i != n) {
      # find index of pivot row (row with max(abs(a[i:n, i]))) using findPivotRowIndex helper function
      pivotRowIndex = findPivotRowIndex(coeffMatrix, i, n)
      
      # if no unique solution
      if(coeffMatrix[pivotRowIndex, i] == 0) return(NA)
      
      # partial pivoting, using swap helper function if current pivot row is not the correct pivot row
      if(i != pivotRowIndex) coeffMatrix = swap(coeffMatrix, i, pivotRowIndex)
    }
    coeffMatrix[i, ] = coeffMatrix[i, ] / coeffMatrix[i, i]     # normalize pivot row (pivot row/pivot element)
    for(j in 1:n) {   
      if(i == j) next
      normalizedRow = coeffMatrix[j, i] * coeffMatrix[i, ]  # get normalized row (value to elim * pivot row)
      coeffMatrix[j, ] = coeffMatrix[j, ] - normalizedRow  # original row - normalized row
    }
    # check main diagonal after normalizing column
    if(!coeffMatrix[i, i] %in% 1) return(NA)      # if main diagonal is not == 1, use %in% to handle NaN values
  }
  
  # build solution vector
  solution = c() # create empty solutions vector
  for(i in 1:n) {
    solution = c(solution, coeffMatrix[i, n+1])
  }
  
  return(list(augcoeffmatrix = coeffMatrix, solution = solution))
}


# ----- Main Code for DPS() ----- #

DPS <- function(df, selectedFood) {
  n = length(selectedFood)
  
  # Step 1: Matrix Construction
  matrix = constructMatrix(df, selectedFood, n)
  
  #print("----- Matrix Construction -----")
  #print(matrix)
  
  # Step 2: Transpose of the Matrix & DUal Problem
  transposedMatrix = t(matrix)
  
  #print("----- Transposed Matrix -----")
  #print(transposedMatrix)
  
  # Step 4: Set-up Initial Tableau
  initialTableau = setupInitialTableau(transposedMatrix, n)
  
  #print("----- Initial Tableau -----")
  #print(initialTableau)
  
  # Step 4: Solving
  solutionsList = solveMaximization(initialTableau)
  if(!is.list(solutionsList)) {
    return(NA)
  }
  
  
  # Generate final output table
  foodServings = c()
  for(i in 1:n) {
    foodServings = c(foodServings, paste("x", i, sep=""))
  }
  foodServings = solutionsList$finalSoln[foodServings]
  foodCosts = df[selectedFood,2]*foodServings
  finalOutputTable = data.frame(Food = df[selectedFood,1],
                                Servings = round(foodServings,2),
                                "Cost($s)"= round(foodCosts,2),
                                row.names = NULL)
  finalOutputTable = subset(finalOutputTable, Servings != 0)
  colnames(finalOutputTable) = c("Food", "Servings", "Costs($)")    # fix column names for special characters i.e. ($)
  
  solutionsList = c(solutionsList, list(finalOutputTable = finalOutputTable, numFood = n))
  
  return(solutionsList)
}


constructMatrix <- function(df, selectedFood, n) {
  
  # Step 1: Matrix Construction
  # Add constraints to matrix
  minNutrientConstraints = c(2000,0,0,0,0,25,50,5000,50,800,10)
  maxNutrientConstraints = c(2250,300,65,2400,300,100,100,50000,20000,1600,30)
  
  matrix = matrix(NA, nrow=23+n, ncol=n+1, byrow=TRUE)
  nutrients = c()
  rowCounter = 1
  for(j in 4:14) {
    # for each nutrient add both min and max constraints to matrix
    for(i in 1:n) {     # iterate over each type of nutrient of food and add it to nutrient, i = nutrient value of food, j = type of nutrient (calories, cholesterol, etc.)
      nutrients = c(nutrients, df[selectedFood[i],j])
    }
    matrix[rowCounter,] = c(nutrients, minNutrientConstraints[j-3])
    matrix[rowCounter+1,] = c(-nutrients, -maxNutrientConstraints[j-3])
    rowCounter = rowCounter + 2
    
    nutrients = c()
  }
  
  for(i in rowCounter:(nrow(matrix)-1)) {
    for(j in 1:ncol(matrix)) {
      if(i-22 == j) matrix[i, j] = -1
      else if(j == ncol(matrix)) matrix[i, j] = -10
      else matrix[i, j] = 0
    }
  }
  
  # Add objective function
  matrix[nrow(matrix),] = c(df[selectedFood,2], 1)    # add unit prices (2nd column of df) of selected food and 1 for the Z at the last row of the matrix
  
  return(matrix)
}


setupInitialTableau <- function(transposedMatrix, n) {
  numCols = ncol(transposedMatrix)-1
  colNames = c()
  for(i in 1:(22+n)) {
    colNames = c(colNames, paste("S", i, sep=""))
  }
  for(i in 1:n) {
    colNames = c(colNames, paste("x", i, sep=""))
  }
  rowNames = c(paste(1:nrow(transposedMatrix)))
  colNames = c(colNames, "Z", "Solution")
  # Step 4: Initial Tableau
  initialTableau = matrix(0, nrow=nrow(transposedMatrix), ncol=numCols+n+2, dimnames=list(rowNames, colNames))
  for(row in 1:(n+1)) {
    for(col in 1:numCols) {
      if(row == n+1) {
        initialTableau[row,col] = -transposedMatrix[row,col]
      } else {
        initialTableau[row,col] = transposedMatrix[row,col]
      }
    }
    initialTableau[row,row+numCols] = 1
    initialTableau[row,ncol(initialTableau)] = transposedMatrix[row,ncol(transposedMatrix)]
  }
  initialTableau[nrow(initialTableau), ncol(initialTableau)] = 0
  return(initialTableau)
}


solveMaximization <- function(initialTableau) {
  initTableau = initialTableau
  n = nrow(initialTableau)
  PC = findHighNegMag(initialTableau[nrow(initialTableau),])
  tableauList = list()
  
  # for basic solution
  basicSoln = list()
  # column names for basic solution
  colNames = c()
  for(i in 1:(22+n-1)) {
    colNames = c(colNames, paste("S", i, sep=""))
  }
  for(i in 1:(n-1)) {
    colNames = c(colNames, paste("x", i, sep=""))
  }
  colNames = c(colNames, "Z")
  
  while(PC != -1) {
    PR = findPR(initialTableau[,PC], initialTableau[,ncol(initialTableau)])
    if(PR == -1) return(NA)
    if(initialTableau[PR, PC] == 0) return(-1)
    initialTableau[PR, ] = initialTableau[PR, ] / initialTableau[PR, PC]     # normalize pivot row (pivot row/pivot element)
    for(j in 1:n) {   
      if(PR == j) next
      normalizedRow = initialTableau[j, PC] * initialTableau[PR, ]  # get normalized row (value to elim * pivot row)
      initialTableau[j, ] = initialTableau[j, ] - normalizedRow  # original row - normalized row
    }
    tableauList[[length(tableauList)+1]] = initialTableau
    PC = findHighNegMag(initialTableau[nrow(initialTableau),])
    
    # for basic solution
    iterSoln = matrix(NA, nrow=1, ncol=ncol(initialTableau)-1, dimnames=list(c(1),colNames[1:length(colNames)]))
    for(col in 1:(ncol(initialTableau)-1)) {
      rowIndex = -1
      for(row in 1:nrow(initialTableau)) {
        if(initialTableau[row, col] == 1) rowIndex = row
        else if(!isTRUE(all.equal(initialTableau[row, col], 0))) rowIndex = -1
      }
      if(rowIndex == -1) iterSoln[1,col] = 0
      else iterSoln[1,col] = initialTableau[rowIndex, ncol(initialTableau)]
    }
    basicSoln[[length(basicSoln)+1]] = iterSoln
  }
  
  # replace basic solution of last iteration, take only the last row and replace z with value at last row, last column
  basicSoln[[length(basicSoln)]] = initialTableau[nrow(initialTableau), 1:(ncol(initialTableau)-1)]
  zVal = initialTableau[nrow(initialTableau), ncol(initialTableau)]
  basicSoln[[length(basicSoln)]]["Z"] = zVal
  
  finalSoln = basicSoln[[length(basicSoln)]]
  
  return(list(initialTableau = initTableau,
              tableaus = tableauList, 
              basicSolns = basicSoln, 
              finalSoln = finalSoln, 
              zVal = zVal))
}


findHighNegMag <- function(bottomRow) {
  negativeIndices = which(bottomRow < 0)    # get indices of negative numbers of bottom row
  negativeBottomRow = c()     # create empty vector
  
  for(i in negativeIndices) {     # iterate over the negative indices
    negativeBottomRow = c(negativeBottomRow, bottomRow[i])    # add the equivalent values to the vector
  }
  
  if(is.null(negativeBottomRow)) return(-1)
  
  highestMagnitude = max(abs(negativeBottomRow))
  PC = which(bottomRow == -highestMagnitude)  # get the highest magnitude among negative
  
  return(PC)    
}


findPR <- function(PC, Sol) {
  n = length(PC)-1
  smallestTR = Inf
  PR = -1
  for(i in 1:n) {
    testRatio = Sol[i] / PC[i]
    if(testRatio > 0 && testRatio < smallestTR) {
      smallestTR = testRatio
      PR = i
    }
  }
  return(PR)
}


# ----- server-side code -----

library(shiny)
library(shinydashboard)
library(DT)

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # ----- QSI code -----
  QSIdata <- eventReactive(input$QSIcalculate,{
    req(input$QSIfile)
    df = read.csv(input$QSIfile$datapath, header=F, col.names=c("x", "f(x)"), check.names = F)
    df = df[order(df[,1]),]
  })
  
  output$QSItable <- renderDataTable({
    req(input$QSIfile)
    datatable(read.csv(input$QSIfile$datapath, header=F, col.names=c("x", "f(x)"), check.names = F))
  })
    
  QSImodel <- eventReactive(input$QSIcalculate,{
    req(input$QSIcalculate)
    QSI(QSIdata(), input$QSIx_estimate)
  })
  
  output$QSIintro_text <- renderUI({
    req(input$QSIcalculate)
    HTML(paste("Your model:", "The f(x) per interval is given by", sep='<br/>'))
  })
  
  output$QSIintervals_output <- renderUI({
    req(input$QSIcalculate)
    intervals = c()
    for(i in 1:length(QSImodel()$intervals)) {
      intervals = c(intervals, paste("[",i,"]", QSImodel()$intervals[i]), "<br/>")
    }
    HTML(intervals)
    #HTML(paste(QSImodel()$intervals, collapse='<br/>'))
  })
  
  output$QSIestimate_output <- renderText({
    req(input$QSIcalculate)
    if(is.numeric(QSImodel()$estimate)) paste("Estimated f(x) = ", round(QSImodel()$estimate, digits=input$QSIround))
    else paste("Estimated f(x) = ", QSImodel()$estimate)
  })
  
  output$QSIfunction_output_text <- renderText({
    req(input$QSIcalculate)
    if(is.numeric(QSImodel()$estimate)) paste("The output function is given by:")
    else paste("")
  })
  
  output$QSIfunction_output <- renderText({
    req(input$QSIcalculate)
    if(is.numeric(QSImodel()$estimate)) paste(QSImodel()$fxn_string)
    else paste("")
  })
  
  output$QSIplot <- renderPlot({
    req(input$QSIcalculate)
    plot(QSIdata(), main="Quadratic Spline Interpolation", pch=10, xlab="x", ylab="f(x)")
    for(i in 1:length(QSImodel()$intervals)) {
      curve(eval(parse(text=QSImodel()$intervals[i]))(x), from=QSIdata()[i,1],to=QSIdata()[i+1,1],col="blue", add=TRUE)
    }
  })
  
  
  # ----- PR code -----
  PRdata <- reactive({
    req(input$PRfile)
    read.csv(input$PRfile$datapath, header=F, col.names=c("x", "f(x)"), check.names = F)
  })
  
  output$PRtable <- renderDataTable({
    req(input$PRfile)
    datatable(PRdata())
  })
  
  PRmodel <- eventReactive(input$PRcalculate,{
    req(input$PRcalculate)
    PolynomialRegression(input$PRdegree, PRdata(), input$PRx_estimate)
  })
  
  isPRSolvable <- eventReactive(input$PRcalculate,{
    if(!is.list(PRmodel())) return(FALSE)
    else return(TRUE)
  })
  
  output$PRplot <- renderPlot({
    req(input$PRcalculate)
    if(isPRSolvable()) {
      plot(PRdata(), main="Polynomial Regression", pch=10, xlab="x", ylab="f(x)")
      curve(PRmodel()$function1(x), col="blue", add=TRUE)
    }
  })
  
  output$PRfunction_output <- renderUI({
    req(input$PRcalculate)
    if(isPRSolvable()) {
      str1 = paste("Your model:")
      str2 = paste("The function f(x) is given by")
      HTML(paste(str1, str2, PRmodel()$fxn_string, sep='<br/>'))
    } else {
      paste("Degree can not be greater than or equal the number of data points!") 
    }
  })
  
  output$PRestimate_output <- renderText({
    req(input$PRcalculate)
    if(isPRSolvable()) paste("Estimated f(x) = ", round(PRmodel()$estimate, digits=input$PRround))
  })
  
  
  # ----- DPS code -----
  
  foodData <- function(){
    read.csv('nutritional_values.csv')
  }
  
  DPSdata <- eventReactive(input$DPScalculate,{
    req(input$foodInput)
    DPSdf = foodData()
    DPSdf[,2] = as.double(substring(DPSdf[,2], 2))
    selectedFoodIndex = c()
    foodNames = foodData()["Foods"]
    for(i in 1:length(input$foodInput)) {
      selectedFoodIndex = c(selectedFoodIndex, which(foodNames == input$foodInput[i]))
    }
    DPS(DPSdf, selectedFoodIndex)
  })
  
  isDPFeasible <- eventReactive(input$DPScalculate,{
    if(!is.list(DPSdata())) return(FALSE)
    else return(TRUE)
  })
  
  output$FoodTable <- renderDataTable({
    datatable(foodData(), options=list(scrollX = TRUE))
  })
  
  output$DPSzVal <- renderText({
    req(input$DPScalculate)
    if(isDPFeasible()) paste("The cost of this optimal diet is $", round(DPSdata()$zVal,2), " per day.", sep="")
    else paste("The problem is infeasible")
  })
  
  output$DPStable <- renderDataTable({
    req(input$DPScalculate)
    if(isDPFeasible()) datatable(DPSdata()$finalOutputTable, options=list(scrollX = TRUE))
  })
  
  output$panels <- renderUI({
    req(input$DPScalculate)
    if(isDPFeasible()) {
      # generate column/row names for basic/final solution
      n = DPSdata()$numFood
      rNames = c(1)
      cNames = c()
      for(i in 1:(22+n)) {
        cNames = c(cNames, paste("S", i, sep=""))
      }
      for(i in 1:n) {
        cNames = c(cNames, paste("x", i, sep=""))
      }
      cNames = c(cNames, "Z")
      
      iterations = length(DPSdata()$tableaus)
      initTableauPanel = list(tabPanel("Initial Tableau", datatable(DPSdata()$initialTableau), options=list(scrollX = TRUE)))
      tabNames = sapply(1:iterations, function(i) c(paste("Iteration", i)))
      captionNames = c()
      captionNames[1:iterations-1] = "Basic Solution"
      captionNames[iterations] = "Final Solution"
      pan = lapply(1:iterations, function(i) { tabPanel(tabNames[i], 
        datatable(DPSdata()$tableaus[[i]], caption=paste("Tableau", i), options=list(scrollX = TRUE)), 
        datatable( matrix(DPSdata()$basicSolns[[i]], ncol=length(DPSdata()$basicSolns[[i]]), dimnames=list(rNames, cNames)) , caption=captionNames[i], options=list(scrollX=TRUE, lengthChange=FALSE, ordering=FALSE, searching=FALSE, paging=FALSE, info=FALSE)) )})
      do.call(tabBox, c(width=12, selected=tabNames[length(tabNames)], initTableauPanel, pan))
    }
  })

})