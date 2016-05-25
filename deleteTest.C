void deleteTest(Int_t x, Int_t y){
  //This is a double matrix in which each element can be a double* (hence an array of arbitrary lenght)
  Double_t ***matrix=new Double_t**[x];

  //Allocation
  for(Int_t i=0; i < x; i++){
    //Here the allocation happens via one "*" less and [] operator
    matrix[i]=new Double_t*[y];
    for (Int_t j = 0; j < y; j++) {
      //Here the allocation happens just allocaing a single double* pointer
      matrix[i][j]=new Double_t;
    }
  }

  //Deallocation must happen bottom top (reverse order)
  for(Int_t i=0; i < x; i++){
    for (Int_t j = 0; j < y; j++) {
      //The inner element has been allocatend simply via "new" so delete is ok
      delete matrix[i][j];
    }
    //The "y" element has been allocated via "[]" so delete [] has to be called ON EACH element separately
    delete[] matrix[i];
  }
  //The last step is to call delete[] on the main pointe ONCE it has been internally unallocated
  delete[] matrix;
  //Please note that delete[][] has no sense: the number of square parenthesis is not related in any way to the depth of the matrix
}
