package TumorRPT;

import java.util.ArrayList;


public class MyUtils {

    public static double[] lastElementOfDoubleArrayList(ArrayList<double[]> myArrayList) {
        if (myArrayList == null || myArrayList.isEmpty()) {
            return null; // or throw an exception based on your requirements
        }
        return myArrayList.get(myArrayList.size() - 1);
    }

	public static boolean isElementPresent(int[] arr, int key) {
		for (int element : arr) {
			if (element == key) {
				return true;
			}
		}
		return false;
	}
}