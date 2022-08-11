//-------------------------------------------------
// Top scope ESI variables
//-------------------------------------------------
int ESIseverityMax {
  value = 0;
  IOstatus = "output";
  description = "
Variable to store maximum severity of the errHandler.ESIs array before
the errHandler is cleared so that it is available to external users.

Example:
  // Could use this in a runfile
  if( ESISeverityMax <  8 ) { cout << \"The point is good\" << endl; }
  if( ESISeverityMax = 8 ) { cout << \"The point is bad\" << endl; }
  if( ESISeverityMax >  8 ) { cout << \"The point is fatal\" << endl; }
";
}

//-------------------------------------------------
//  Function to set all top scope ESI variables
//-------------------------------------------------
void ESIpromote() {
	ESIseverityMax = getESIseverityMax();
}
ESIpromote.description = "
Function sets top scope ESI variables with information from the errHandler.ESIs array.
";

//-------------------------------------------------
// Unit functions for getting individual values
//-------------------------------------------------
int getESIseverityMax() {

  int i,severityMax,severityTmp,ESITemp[];
  
  severityMax = 0;
  ESITemp = .errHandler.ESIs;
  
  for(i=0;i<ESITemp.entries();i++) {
    severityTmp = getESIseverity(ESITemp[i]);
    severityMax = max(severityMax,severityTmp);
  }

  return severityMax;
}
getESIseverityMax.description = "
Function to return maximum severity of the errHandler.ESIs array.
This function is for internal use before the errHandler is cleared.

Example:
  // Could use this inside the solverSequence
  if( getESIseverityMax() <  8 ) { cout << \"The point is good\" << endl; }
  if( getESIseverityMax() = 8 ) { cout << \"The point is bad\" << endl; }
  if( getESIseverityMax() >  8 ) { cout << \"The point is fatal\" << endl; }
";

//-------------------------------------------------
// Constants from the ARP 5571 standard
//-------------------------------------------------
int ESImaxDigits {
	value = 9;
	description = "Maximum number of digits in the ESI code";
	IOstatus = "const";
}
// Note: Additional information on the ARP 5571 ESI format is encoded
//       in the utility parsing functions below.

//-------------------------------------------------
// Utility functions for parsing ESI codes
//-------------------------------------------------
// ARP 5571 standard for ESI codes
// ESI = (S)SPPPCCOO
// Digits   1,2: Severity
// Digits 3,4,5: Purpose
// Digits   6,7: Category
// Digits   8,9: Owner

int getESIseverity(int ESI) {
	return getDigits(ESI,1,2,ESImaxDigits);
}
getESIseverity.description = "Function to return ESI severity code.";

int getESIpurpose(int ESI) {
	return getDigits(ESI,3,5,ESImaxDigits);
}
getESIpurpose.description = "Function to return ESI purpose code.";

int getESIcategory(int ESI) {
	return getDigits(ESI,6,7,ESImaxDigits);
}
getESIcategory.description = "Function to return ESI category code.";

int getESIowner(int ESI) {
	return getDigits(ESI,8,9,ESImaxDigits);
}
getESIowner.description = "Function to return ESI owner code.";

// getDigits() Arguments
//
//        var: integer to be parsed
//  leftDigit: position of the leftmost digit that will be extracted 
//             (1 is the leftmost digit; maxDigits is the rightmost)
// rightDigit: position of the rightmost digit that will be extracted 
//             (1 is the leftmost digit; maxDigits is the rightmost)
//  maxDigits: maximum number of digits in the integer 
//             (integer can be any fewer number of digits, but position
//              inputs are always based on the maximum value).
//
// getDigits() Examples
//
// getDigits(12345,2,4,5) = 234
// getDigits(2345,2,4,5) = 234
//  
int getDigits(int var, int leftDigit, int rightDigit, int maxDigits) {

	// Initialize output assuming invalid input
	int out = NaN;

	// Validate input
	int id = 0002200;
	ESOregCreate(id, 7, getName() + " input error");
		
	if ((floor(log10(var)) + 1) > maxDigits) {
		ESOreport(id,"integer length exceeds maxDigits", TRUE);
	} else if (leftDigit > rightDigit) {
		ESOreport(id,"leftDigit exceeds rightDigit", TRUE);
	} else if (rightDigit > maxDigits) {
		ESOreport(id,"rightDigit exceeds maxDigits", TRUE);
	} else if (leftDigit < 1) {
		ESOreport(id,"leftDigit must be positive", TRUE);
	} else if (var > (2**31 - 2)) {
		ESOreport(id,"integer exceeds maximum size", TRUE);
	} else if (var < 0) {
		ESOreport(id,"integer must be non-negative", TRUE);
	} else {
		// Correct output if input is valid
		// Strip trailing digits
		out = (var/10**(maxDigits-rightDigit));
		// Strip leading digits
		out = out % 10**(rightDigit - leftDigit + 1);
	}
	
	return out;
}
getDigits.description = "Function to return the integer digits located between specified left and right digits (inclusive)";