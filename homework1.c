// C program for Sorting Algorithms
#include <stdio.h>
#include <math.h>
#include <time.h>

/* Constant Definition */
#define ARRAY_SIZE 10000
#define ITERATION_COUNT 10

/* Function to sort an array using insertion sort*/
void insertionSort(int arr[], int n)
{
   int i, key, j;
   for (i = 1; i < n; i++)
   {
       key = arr[i];
       j = i-1;
 
       /* Move elements of arr[0..i-1], that are
          greater than key, to one position ahead
          of their current position */
       while (j >= 0 && arr[j] > key)
       {
           arr[j+1] = arr[j];
           j = j-1;
       }
       arr[j+1] = key;
   }
}


// Swap Function
void swap(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}
 
// Selection Sort Function
void selectionSort(int arr[], int n)
{
    int i, j, min_idx;
 
    // One by one move boundary of unsorted subarray
    for (i = 0; i < n-1; i++)
    {
        // Find the minimum element in unsorted array
        min_idx = i;
        for (j = i+1; j < n; j++)
          if (arr[j] < arr[min_idx])
            min_idx = j;
 
        // Swap the found minimum element with the first element
        swap(&arr[min_idx], &arr[i]);
    }
}
 

void bubbleSort(int arr[], int n)
{
   int i,j;	 
   int sorted = 0; 
   int last = n-1;
      
   for (i = 0; (i < last) && (sorted==0); i++){
      sorted = 1;
      for (j=last; j > i; j--)
         if (arr[j-1] > arr[j]){
            swap(&arr[j],&arr[j-1]);
            sorted = 0; // signal exchange
         }
    }
}  

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(int arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 =  r - m;
 
    /* create temp arrays */
    int L[n1], R[n2];
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

// Merge Sort Function 
/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort(int arr[], int l, int r)
{
    if (l < r)
    {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l+(r-l)/2;
 
        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m+1, r);
 
        merge(arr, l, m, r);
    }
}


//quick Sort function to Sort Integer array list
void quickSort(int array[], int firstIndex, int lastIndex)
{
    //declaaring index variables
    int pivotIndex, temp, index1, index2;

    if(firstIndex < lastIndex)
    {
        //assigninh first element index as pivot element
        pivotIndex = firstIndex;
        index1 = firstIndex;
        index2 = lastIndex;

        //Sorting in Ascending order with quick sort
        while(index1 < index2)
        {
            while(array[index1] <= array[pivotIndex] && index1 < lastIndex)
            {
                index1++;
            }
            while(array[index2]>array[pivotIndex])
            {
                index2--;
            }

            if(index1<index2)
            {
                //Swapping opertation
                temp = array[index1];
                array[index1] = array[index2];
                array[index2] = temp;
            }
        }

        //At the end of first iteration, swap pivot element with index2 element
        temp = array[pivotIndex];
        array[pivotIndex] = array[index2];
        array[index2] = temp;

        //Recursive call for quick sort, with partiontioning
        quickSort(array, firstIndex, index2-1);
        quickSort(array, index2+1, lastIndex);
    }
}

// To heapify a subtree rooted with node i which is
// an index in arr[]. n is size of heap
void heapify(int arr[], int n, int i)
{
    int smallest = i;  // Initialize smallest as root
    int l = 2*i + 1;  // left = 2*i + 1
    int r = 2*i + 2;  // right = 2*i + 2
 
    // If left child is larger than root
    if (l < n && arr[l] < arr[smallest])
        smallest = l;
 
    // If right child is larger than largest so far
    if (r < n && arr[r] < arr[smallest])
        smallest = r;
 
    // If largest is not root
    if (smallest != i)
    {
        swap(&arr[i], &arr[smallest]);
 
        // Recursively heapify the affected sub-tree
        heapify(arr, n, smallest);
    }
}


// main function to do heap sort
void heapSort(int arr[], int n)
{
    // Build heap (rearrange array)
    int i; 
    for (i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);
 
    // One by one extract an element from heap
    for (i=n-1; i>=0; i--)
    {
        // Move current root to end
        swap(&arr[0], &arr[i]);
 
        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    }
}

  
// A utility function ot print an array of size n
void printArray(int arr[], int n)
{
   int i;
   for (i=0; i < n; i++)
       printf("%d ", arr[i]);
   printf("\n");
}
 
 
// Three array as Sorted Order, Reverse Order and Random Order
// Original Array
int sorted_order[ARRAY_SIZE];
int reverse_order[ARRAY_SIZE];
int random_order[ARRAY_SIZE];

// Enumeration types for stroing passed time values
enum sort_types   {Insertion,Selection,Bubble,Merge,Quick,Heap};
enum sort_cases   {Best,Worst,Average};
enum min_max_time {Min,Max};

double average_time_result[6][3];
double minmax_time_result[6][3][2];

// This will be used temporary array
int arr[ARRAY_SIZE];


void prepare_array()
{
	int i=0;
	int j=0;
	int k=0;
	
	// Prepare sorted_order array
	for(i=0;i<ARRAY_SIZE;i++)
	{
		sorted_order[i]=i+1;
	}
	
    // Prepare reverse_order array
	for(i=0;i<ARRAY_SIZE;i++)
	{
		reverse_order[i]=ARRAY_SIZE-i;
	}
	
    // Prepare random_order array
    srand(time(NULL)); //required for "randomness"
	for(i=0;i<ARRAY_SIZE;i++)
	{
		random_order[i]=rand()%ARRAY_SIZE;
	}
	
	// Initialize temporary array
	for(i=0;i<ARRAY_SIZE;i++)
	{
		arr[i]=0;
	}
	
	// Initialize average_time_result matrix
	for(i=0;i<6;i++)
	{
		for(j=0;j<3;j++)
		{
			average_time_result[i][j]=0;
			minmax_time_result[i][j][Min]=ARRAY_SIZE;
	        minmax_time_result[i][j][Max]=0;
		}
		
	}
	
}
	 
void copy_array(int source[ARRAY_SIZE])
{
	int i=0; 
	for(i=0;i<ARRAY_SIZE;i++)
    {
    	arr[i]=source[i];
	}
}	 


void test_insertion_best_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++) // Algorithm will be run for ITERATION count time have safe result
	{
		copy_array(sorted_order); // Fill arr[] with sorted_order
    	t = clock();
    	insertionSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;  // Time consumed during sort
    	average_time_result[Insertion][Best]=average_time_result[Insertion][Best] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Insertion][Best][Min])
    	{
    		minmax_time_result[Insertion][Best][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Insertion][Best][Max])
    	{
    		minmax_time_result[Insertion][Best][Max]=pased_time;
		}
    }

    printf("Insertion Sort Best Case Average Time: %f  \n", average_time_result[Insertion][Best]/ITERATION_COUNT);
    printf("Insertion Sort Best Case Minimum Time: %f  \n", minmax_time_result[Insertion][Best][Min]);
    printf("Insertion Sort Best Case Maximum Time: %f  \n", minmax_time_result[Insertion][Best][Max]);

    printf("\n"); 

}	 

void test_insertion_worst_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(reverse_order);
    	t = clock();
    	insertionSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Insertion][Worst]=average_time_result[Insertion][Worst] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Insertion][Worst][Min])
    	{
    		minmax_time_result[Insertion][Worst][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Insertion][Worst][Max])
    	{
    		minmax_time_result[Insertion][Worst][Max]=pased_time;
		}
    }

    printf("Insertion Sort Worst Case Average Time: %f  \n", average_time_result[Insertion][Worst]/ITERATION_COUNT);
    printf("Insertion Sort Worst Case Minimum Time: %f  \n", minmax_time_result[Insertion][Worst][Min]);
    printf("Insertion Sort Worst Case Maximum Time: %f  \n", minmax_time_result[Insertion][Worst][Max]);

    printf("\n"); 
    
}

void test_insertion_average_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(random_order);
    	t = clock();
    	insertionSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Insertion][Average]=average_time_result[Insertion][Average] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Insertion][Average][Min])
    	{
    		minmax_time_result[Insertion][Average][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Insertion][Average][Max])
    	{
    		minmax_time_result[Insertion][Average][Max]=pased_time;
		}
		
    }
    
    printf("Insertion Sort Average Case Average Time: %f  \n", average_time_result[Insertion][Average]/ITERATION_COUNT);
    printf("Insertion Sort Average Case Minimum Time: %f  \n", minmax_time_result[Insertion][Average][Min]);
    printf("Insertion Sort Average Case Maximum Time: %f  \n", minmax_time_result[Insertion][Average][Max]);

    printf("\n"); 


}

void all_insertion_sort_cases()
{
	test_insertion_best_case();
	test_insertion_worst_case();
	test_insertion_average_case();
}

void test_selection_best_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(sorted_order);
    	t = clock();
    	selectionSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Selection][Best]=average_time_result[Selection][Best] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Selection][Best][Min])
    	{
    		minmax_time_result[Selection][Best][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Selection][Best][Max])
    	{
    		minmax_time_result[Selection][Best][Max]=pased_time;
		}
    }

    printf("Selection Sort Best Case Average Time: %f  \n", average_time_result[Selection][Best]/ITERATION_COUNT);
    printf("Selection Sort Best Case Minimum Time: %f  \n", minmax_time_result[Selection][Best][Min]);
    printf("Selection Sort Best Case Maximum Time: %f  \n", minmax_time_result[Selection][Best][Max]);

    printf("\n"); 


}	 

void test_selection_worst_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(reverse_order);
    	t = clock();
    	selectionSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Selection][Worst]=average_time_result[Selection][Worst] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Selection][Worst][Min])
    	{
    		minmax_time_result[Selection][Worst][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Selection][Worst][Max])
    	{
    		minmax_time_result[Selection][Worst][Max]=pased_time;
		}
    }
    
    printf("Selection Sort Worst Case Average Time: %f  \n", average_time_result[Selection][Worst]/ITERATION_COUNT);
    printf("Selection Sort Worst Case Minimum Time: %f  \n", minmax_time_result[Selection][Worst][Min]);
    printf("Selection Sort Worst Case Maximum Time: %f  \n", minmax_time_result[Selection][Worst][Max]);

    printf("\n"); 



}

void test_selection_average_case()
{
	int i=0;
	clock_t t; // Time variable
    double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(random_order);
    	t = clock();
    	selectionSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Selection][Average]=average_time_result[Selection][Average] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Selection][Average][Min])
    	{
    		minmax_time_result[Selection][Average][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Selection][Average][Max])
    	{
    		minmax_time_result[Selection][Average][Max]=pased_time;
		}
    }
    
    printf("Selection Sort Average Case Average Time: %f  \n", average_time_result[Selection][Average]/ITERATION_COUNT);
    printf("Selection Sort Average Case Minimum Time: %f  \n", minmax_time_result[Selection][Average][Min]);
    printf("Selection Sort Average Case Maximum Time: %f  \n", minmax_time_result[Selection][Average][Max]);

    printf("\n"); 


}

void all_selection_sort_cases()
{
	test_selection_best_case();
	test_selection_worst_case();
	test_selection_average_case();
}

void test_bubble_best_case()
{
	int i=0;
	clock_t t; // Time variable
    double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(sorted_order);
    	t = clock();
    	bubbleSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Bubble][Best]=average_time_result[Bubble][Best] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Bubble][Best][Min])
    	{
    		minmax_time_result[Bubble][Best][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Bubble][Best][Max])
    	{
    		minmax_time_result[Bubble][Best][Max]=pased_time;
		}
    }
    
    printf("Bubble Sort Best Case Average Time: %f  \n", average_time_result[Bubble][Best]/ITERATION_COUNT);
    printf("Bubble Sort Best Case Minimum Time: %f  \n", minmax_time_result[Bubble][Best][Min]);
    printf("Bubble Sort Best Case Maximum Time: %f  \n", minmax_time_result[Bubble][Best][Max]);

    printf("\n"); 


}	 

void test_bubble_worst_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(reverse_order);
    	t = clock();
    	bubbleSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Bubble][Worst]=average_time_result[Bubble][Worst] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Bubble][Worst][Min])
    	{
    		minmax_time_result[Bubble][Worst][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Bubble][Worst][Max])
    	{
    		minmax_time_result[Bubble][Worst][Max]=pased_time;
		}
    }
    
    printf("Bubble Sort Worst Case Average Time: %f  \n", average_time_result[Bubble][Worst]/ITERATION_COUNT);
    printf("Bubble Sort Worst Case Minimum Time: %f  \n", minmax_time_result[Bubble][Worst][Min]);
    printf("Bubble Sort Worst Case Maximum Time: %f  \n", minmax_time_result[Bubble][Worst][Max]);

    printf("\n"); 
 

}

void test_bubble_average_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(random_order);
    	t = clock();
    	bubbleSort(arr, ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Bubble][Average]=average_time_result[Bubble][Average] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Bubble][Average][Min])
    	{
    		minmax_time_result[Bubble][Average][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Bubble][Average][Max])
    	{
    		minmax_time_result[Bubble][Average][Max]=pased_time;
		}
    }
    
    printf("Bubble Sort Average Case Average Time: %f  \n", average_time_result[Bubble][Average]/ITERATION_COUNT);
    printf("Bubble Sort Average Case Minimum Time: %f  \n", minmax_time_result[Bubble][Average][Min]);
    printf("Bubble Sort Average Case Maximum Time: %f  \n", minmax_time_result[Bubble][Average][Max]);

    printf("\n"); 


}

void all_bubble_sort_cases()
{
	test_bubble_best_case();
	test_bubble_worst_case();
	test_bubble_average_case();
}

void test_merge_best_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(sorted_order);
    	t = clock();
    	mergeSort(arr, 0,ARRAY_SIZE-1);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Merge][Best]=average_time_result[Merge][Best] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Merge][Best][Min])
    	{
    		minmax_time_result[Merge][Best][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Merge][Best][Max])
    	{
    		minmax_time_result[Merge][Best][Max]=pased_time;
		}
    }
    
    printf("Merge Sort Best Case Average Time: %f  \n", average_time_result[Merge][Best]/ITERATION_COUNT);
    printf("Merge Sort Best Case Minimum Time: %f  \n", minmax_time_result[Merge][Best][Min]);
    printf("Merge Sort Best Case Maximum Time: %f  \n", minmax_time_result[Merge][Best][Max]);

    printf("\n"); 



}	 

void test_merge_worst_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(reverse_order);
    	t = clock();
    	mergeSort(arr, 0,ARRAY_SIZE-1);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Merge][Worst]=average_time_result[Merge][Worst] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Merge][Worst][Min])
    	{
    		minmax_time_result[Merge][Worst][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Merge][Worst][Max])
    	{
    		minmax_time_result[Merge][Worst][Max]=pased_time;
		}
    }
    
    printf("Merge Sort Worst Case Average Time: %f  \n", average_time_result[Merge][Worst]/ITERATION_COUNT);
    printf("Merge Sort Worst Case Minimum Time: %f  \n", minmax_time_result[Merge][Worst][Min]);
    printf("Merge Sort Worst Case Maximum Time: %f  \n", minmax_time_result[Merge][Worst][Max]);

    printf("\n"); 
 

}

void test_merge_average_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(random_order);
    	t = clock();
    	mergeSort(arr, 0,ARRAY_SIZE-1);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Merge][Average]=average_time_result[Merge][Average] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Merge][Average][Min])
    	{
    		minmax_time_result[Merge][Average][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Merge][Average][Max])
    	{
    		minmax_time_result[Merge][Average][Max]=pased_time;
		}
    }
    
    printf("Merge Sort Average Case Average Time: %f  \n", average_time_result[Merge][Average]/ITERATION_COUNT);
    printf("Merge Sort Average Case Minimum Time: %f  \n", minmax_time_result[Merge][Average][Min]);
    printf("Merge Sort Average Case Maximum Time: %f  \n", minmax_time_result[Merge][Average][Max]);

    printf("\n"); 

}

void all_merge_sort_cases()
{
	test_merge_best_case();
	test_merge_worst_case();
	test_merge_average_case();
}

void test_quick_best_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(sorted_order);
    	t = clock();
    	quickSort(arr, 0,ARRAY_SIZE-1);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Quick][Best]=average_time_result[Quick][Best] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Quick][Best][Min])
    	{
    		minmax_time_result[Quick][Best][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Quick][Best][Max])
    	{
    		minmax_time_result[Quick][Best][Max]=pased_time;
		}
    }
 
    printf("Quick Sort Best Case Average Time: %f  \n", average_time_result[Quick][Best]/ITERATION_COUNT);
    printf("Quick Sort Best Case Minimum Time: %f  \n", minmax_time_result[Quick][Best][Min]);
    printf("Quick Sort Best Case Maximum Time: %f  \n", minmax_time_result[Quick][Best][Max]);

    printf("\n"); 

}	 

void test_quick_worst_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(reverse_order);
    	t = clock();
    	quickSort(arr, 0,ARRAY_SIZE-1);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Quick][Worst]=average_time_result[Quick][Worst] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Quick][Worst][Min])
    	{
    		minmax_time_result[Quick][Worst][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Quick][Worst][Max])
    	{
    		minmax_time_result[Quick][Worst][Max]=pased_time;
		}
    }
    
    printf("Quick Sort Worst Case Average Time: %f  \n", average_time_result[Quick][Worst]/ITERATION_COUNT);
    printf("Quick Sort Worst Case Minimum Time: %f  \n", minmax_time_result[Quick][Worst][Min]);
    printf("Quick Sort Worst Case Maximum Time: %f  \n", minmax_time_result[Quick][Worst][Max]);

    printf("\n"); 


}

void test_quick_average_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(random_order);
    	t = clock();
    	quickSort(arr, 0,ARRAY_SIZE-1);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Quick][Average]=average_time_result[Quick][Average] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Quick][Average][Min])
    	{
    		minmax_time_result[Quick][Average][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Quick][Average][Max])
    	{
    		minmax_time_result[Quick][Average][Max]=pased_time;
		}
    }

    printf("Quick Sort Average Case Average Time: %f  \n", average_time_result[Quick][Average]/ITERATION_COUNT);
    printf("Quick Sort Average Case Minimum Time: %f  \n", minmax_time_result[Quick][Average][Min]);
    printf("Quick Sort Average Case Maximum Time: %f  \n", minmax_time_result[Quick][Average][Max]);

    printf("\n"); 

}

void all_quick_sort_cases()
{
	test_quick_best_case();
	test_quick_worst_case();
	test_quick_average_case();
}

void test_heap_best_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(sorted_order);
    	t = clock();
    	heapSort(arr,ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Heap][Best]=average_time_result[Heap][Best] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Heap][Best][Min])
    	{
    		minmax_time_result[Heap][Best][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Heap][Best][Max])
    	{
    		minmax_time_result[Heap][Best][Max]=pased_time;
		}
    }

    printf("Heap Sort Best Case Average Time: %f  \n", average_time_result[Heap][Best]/ITERATION_COUNT);
    printf("Heap Sort Best Case Minimum Time: %f  \n", minmax_time_result[Heap][Best][Min]);
    printf("Heap Sort Best Case Maximum Time: %f  \n", minmax_time_result[Heap][Best][Max]);

    printf("\n"); 

}	 

void test_heap_worst_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(reverse_order);
    	t = clock();
    	heapSort(arr,ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Heap][Worst]=average_time_result[Heap][Worst] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Heap][Worst][Min])
    	{
    		minmax_time_result[Heap][Worst][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Heap][Worst][Max])
    	{
    		minmax_time_result[Heap][Worst][Max]=pased_time;
		}
    }
    
    printf("Heap Sort Worst Case Average Time: %f  \n", average_time_result[Heap][Worst]/ITERATION_COUNT);
    printf("Heap Sort Worst Case Minimum Time: %f  \n", minmax_time_result[Heap][Worst][Min]);
    printf("Heap Sort Worst Case Maximum Time: %f  \n", minmax_time_result[Heap][Worst][Max]);

    printf("\n"); 

}

void test_heap_average_case()
{
	int i=0;
	clock_t t; // Time variable
	double pased_time=0;
	for(i=0;i<ITERATION_COUNT;i++)
	{
		copy_array(random_order);
    	t = clock();
    	heapSort(arr,ARRAY_SIZE);
    	t = clock() - t;
    	pased_time=(double)t/CLOCKS_PER_SEC;
    	average_time_result[Heap][Average]=average_time_result[Heap][Average] + pased_time;
    	
    	// Update Minimum Time
    	if (pased_time < minmax_time_result[Heap][Average][Min])
    	{
    		minmax_time_result[Heap][Average][Min]=pased_time;
		}
		
		// Update Maximum Time
		if (pased_time > minmax_time_result[Heap][Average][Max])
    	{
    		minmax_time_result[Heap][Average][Max]=pased_time;
		}
    }

    printf("Heap Sort Average Case Average Time: %f  \n", average_time_result[Heap][Average]/ITERATION_COUNT);
    printf("Heap Sort Average Case Minimum Time: %f  \n", minmax_time_result[Heap][Average][Min]);
    printf("Heap Sort Average Case Maximum Time: %f  \n", minmax_time_result[Heap][Average][Max]);

    printf("\n"); 
    
}

void all_heap_sort_cases()
{
	test_heap_best_case();
	test_heap_worst_case();
	test_heap_average_case();
}


int main()
{
	
    prepare_array(); // Prerequisite work
    
    all_insertion_sort_cases(); // Calling Insertion Sort
    
    all_selection_sort_cases(); // Calling Selection Sort
    
    all_bubble_sort_cases(); // Calling Bubble Sort
    
    all_merge_sort_cases(); // Calling Merge Sort
    
    all_quick_sort_cases(); // Calling Quick Sort
    
    all_heap_sort_cases(); // Calling Heap Sort
    
    return 0;
} 
