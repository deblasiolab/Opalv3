package opal.align;

public class SequenceIdPair {
	
	public static String makeString (int a, int b) {

		
//		if (a<b) {
			return a + "," + b;
//		} else {
//			return b + "," + a;
//		}
	}

	public static int[] splitString (String s) {
		String[] nums = s.split(",");
		
		int[] ret = new int[2];
		
		
		ret[0] =  new Integer(nums[0]).intValue();
		ret[1] =  new Integer(nums[1]).intValue();
		return ret;
	}
	
}
