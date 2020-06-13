package opal.exceptions;

/**
 * Exception class for Opal. If it should crash ... 
 * throw this so the wrapper can catch it 
 * @author Travis Wheeler
 */
public class GenericOpalException extends RuntimeException {
    /**
     * Construct this exception object.
     * @param message the error message.
     */
    public GenericOpalException( String message ) {
        super( message );
    }
}