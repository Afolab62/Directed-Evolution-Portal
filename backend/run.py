from flask import Flask
from flask_cors import CORS
from flask_session import Session
from config import Config
from routes.auth import auth_bp
from routes.experiments import experiments_bp
from routes.uniprot import uniprot_bp
from database import init_db


def create_app():
    """Create and configure the Flask application"""
    app = Flask(__name__)
    app.config.from_object(Config)
    
    # Initialize extensions
    CORS(app, origins=app.config['CORS_ORIGINS'], supports_credentials=True)
    Session(app)
    
    # Initialize database
    with app.app_context():
        init_db()
    
    # Register blueprints
    app.register_blueprint(auth_bp)
    app.register_blueprint(experiments_bp)
    app.register_blueprint(uniprot_bp)
    
    # Health check endpoint
    @app.route('/health', methods=['GET'])
    def health():
        return {'status': 'ok'}, 200
    
    return app


if __name__ == '__main__':
    app = create_app()
    app.run(host='0.0.0.0', port=8000, debug=True)
