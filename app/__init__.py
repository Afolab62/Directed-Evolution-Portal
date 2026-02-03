from flask import Flask

from .config import Config
from .extensions import db, migrate, login_manager


def create_app():
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_object(Config)

    # init extensions
    db.init_app(app)
    migrate.init_app(app, db)
    login_manager.init_app(app)

    # register blueprints
    from .routes.auth import auth_bp
    from .routes.main import main_bp
    from .routes.staging import staging_bp
    app.register_blueprint(auth_bp)
    app.register_blueprint(main_bp)
    app.register_blueprint(staging_bp)

    return app
