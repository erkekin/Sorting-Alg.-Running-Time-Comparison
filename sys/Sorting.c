#include <stdio.h>
#include <stdlib.h>
#include "Sorting.h"
#include <limits.h>
#include <stdlib.h>

typedef unsigned uint;
#define swap(a, b) { tmp = a; a = b; b = tmp; }
#define each(i, x) for (i = 0; i < x; i++)


/****************************************************
 * Uses bubble sort to sort A[0..N-1]
 ****************************************************/
void bubble_sort(int a[],int n)
{
    int i,j,temp;

    for(i=1;i< n;i++)
    {
        for(j=0;j< n-1;j++)
            if(a[j]>a[j+1])
            {
                temp=a[j];
                a[j]=a[j+1];
                a[j+1]=temp;
            }
    }
}
void BubbleSort(int A[], int N){

    bubble_sort(A, N);

    
} //end-BubbleSort

/****************************************************
 * Uses selection sort to sort A[0..N-1]
 ****************************************************/
void selectSort(int arr[], int n)
{
    //pos_min is short for position of min
	int pos_min,temp;
    
	for (int i=0; i < n-1; i++)
	{
	    pos_min = i;//set pos_min to the current index of array
		
		for (int j=i+1; j < n; j++)
		{
            
            if (arr[j] < arr[pos_min])
                pos_min=j;
            //pos_min will keep track of the index that min is in, this is needed when a swap happens
		}
		
        //if pos_min no longer equals i than a smaller value must have been found, so a swap must occur
        if (pos_min != i)
        {
            temp = arr[i];
            arr[i] = arr[pos_min];
            arr[pos_min] = temp;
        }
	}
}
void SelectionSort(int A[], int N){
    selectSort(A, N);
} //end-SelectionSort

/****************************************************
 * Uses insertion sort to sort A[0..N-1]
 ****************************************************/
void InsertionSort1(int *array, int n)
{
    int c, d, t;
 
    
    for (c = 1 ; c <= n - 1; c++) {
        d = c;
        
        while ( d > 0 && array[d] < array[d-1]) {
            t          = array[d];
            array[d]   = array[d-1];
            array[d-1] = t;
            
            d--;
        }
    }

    
}
void InsertionSort(int A[], int N){
    InsertionSort1(A, N);
} //end-InsertionSort

/****************************************************
 * Uses mergesort to sort A[0..N-1]
 ****************************************************/
void Merge1(int *array, int left, int mid, int right)
{
    /*We need a Temporary array to store the new sorted part*/
    int tempArray[right-left+1];
    int pos=0,lpos = left,rpos = mid + 1;
    while(lpos <= mid && rpos <= right)
    {
        if(array[lpos] < array[rpos])
        {
            tempArray[pos++] = array[lpos++];
        }
        else
        {
            tempArray[pos++] = array[rpos++];
        }
    }
    while(lpos <= mid)  tempArray[pos++] = array[lpos++];
    while(rpos <= right)tempArray[pos++] = array[rpos++];
    int iter;
    /* Copy back the sorted array to the original array */
    for(iter = 0;iter < pos; iter++)
    {
        array[iter+left] = tempArray[iter];
    }
    return;
}
void MergeSort1(int *array, int left, int right)
{
    int mid = (left+right)/2;
    /* We have to sort only when left<right because when left=right it is anyhow sorted*/
    if(left<right)
    {
        /* Sort the left part */
        MergeSort1(array,left,mid);
        /* Sort the right part */
        MergeSort1(array,mid+1,right);
        /* Merge the two sorted parts */
        Merge1(array,left,mid,right);
    }
}
/* Merge functions merges the two sorted parts. Sorted parts will be from [left, mid] and [mid+1, right].
 */

void MergeSort(int A[], int N){
    MergeSort1(A, 0, N-1);
} //end-MergeSort

/****************************************************
 * Uses quicksort sort to sort A[0..N-1]
 ****************************************************/
void quicksort(int x[500],int first,int last){
    int pivot,j,temp,i;
    
    if(first<last){
        pivot=first;
        i=first;
        j=last;
        
        while(i<j){
            while(x[i]<=x[pivot]&&i<last)
                i++;
            while(x[j]>x[pivot])
                j--;
            if(i<j){
                temp=x[i];
                x[i]=x[j];
                x[j]=temp;
            }
        }
        
        temp=x[pivot];
        x[pivot]=x[j];
        x[j]=temp;
        quicksort(x,first,j-1);
        quicksort(x,j+1,last);
        
    }
}
void QuickSort(int A[], int N){
    quicksort(A, 0, N-1);
} //end-QuickSort

/****************************************************
 * Uses heapsort to sort A[0..N]
 * NOTICE: The first element is in A[0] not in A[1]
 ****************************************************/



void heapify(int array[], int n)
{
    int item,i,j,k;
    
    for(k=1 ; k<n ; k++)
    {
        item = array[k];
        i = k;
        j = (i-1)/2;
        
        while( (i>0) && (item>array[j]) )
        {
            array[i] = array[j];
            i = j;
            j = (i-1)/2;
        }
        array[i] = item;
    }
}

void adjust(int array[], int n)
{
    int item,i,j;
    
    j = 0;
    item = array[j];
    i = 2*j+1;
    
    while(i<=n-1)
    {
        if(i+1 <= n-1)
            if(array[i] < array[i+1])
                i++;
        if(item < array[i])
        {
            array[j] = array[i];
            j = i;
            i = 2*j+1;
        }
        else
            break;
    }
    array[j] = item;
}
void heapsorting(int array[], int n)
{
    int i,t;
    
    heapify(array,n);
    
    for(i=n-1 ; i>0 ; i--)
    {
        t = array[0];
        array[0] = array[i];
        array[i] = t;
        adjust(array,i);
    }
}
void HeapSort(int A[], int N){
    heapsorting(A, N);
} //end-HeapSort

/****************************************************
 * Uses radixsort to sort A[0..N]
 ****************************************************/

typedef unsigned uint;
#define swap(a, b) { tmp = a; a = b; b = tmp; }
#define each(i, x) for (i = 0; i < x; i++)

/* sort unsigned ints */

static void rad_sort_u(uint *from, uint *to, uint bit)
{
	if (!bit || to < from + 1) return;
    
	uint *ll = from, *rr = to - 1, tmp;
	while (1) {
		/* find left most with bit, and right most without bit, swap */
		while (ll < rr && !(*ll & bit)) ll++;
		while (ll < rr &&  (*rr & bit)) rr--;
		if (ll >= rr) break;
		swap(*ll, *rr);
	}
    
	if (!(bit & *ll) && ll < to) ll++;
	bit >>= 1;
    
	rad_sort_u(from, ll, bit);
	rad_sort_u(ll, to, bit);
}

/* sort signed ints: flip highest bit, sort as unsigned, flip back */
static void radix_sort(int *a, const size_t len)
{
	size_t i;
	uint *x = (uint*) a;
    
	each(i, len) x[i] ^= INT_MIN;
	rad_sort_u(x, x + len, INT_MIN);
	each(i, len) x[i] ^= INT_MIN;
}


void RadixSort(int A[], int N){

    radix_sort(A, N);
    
} //end-RadixSort

/****************************************************
 * Uses countingsort to sort A[0..N]
 ****************************************************/
void counting_sort(int *array, int n)
{
    int i, min, max;
    
    min = max = array[0];
    for(i=1; i < n; i++) {
        if ( array[i] < min ) {
            min = array[i];
        } else if ( array[i] > max ) {
            max = array[i];
        }
    }
}
void counting_sort_mm(int *array, int n, int min, int max)
{
    int i, j, z;
    
    int range = max - min + 1;
    int *count = malloc(range * sizeof(*array));
    
    for(i = 0; i < range; i++) count[i] = 0;
    for(i = 0; i < n; i++) count[ array[i] - min ]++;
    
    for(i = min, z = 0; i <= max; i++) {
        for(j = 0; j < count[i - min]; j++) {
            array[z++] = i;
        }
    } 
    
    free(count);
}

void CountingSort(int A[], int N){

    counting_sort_mm(A, N, 0, 4444444);
  // Fill this in
} //end-CountingSort

